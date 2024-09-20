# Title: Dissertation_Ch3_Data_Analsysis
# Author: Greg Goodrum
# Last update: 8/22/2023
# Contact: greg.goodrum@usu.edu
# Description: Code to produce results used in MANUSCRIPT TITLE HERE

# NOTES:
# Calculating connectivty indices
# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# References:
# Gradient boosting models
# https://www.datatechnotes.com/2018/03/classification-with-gradient-boosting.html
# https://amunategui.github.io/binary-outcome-modeling/
# https://bookdown.org/rehk/stm1001_dsm_introduction_to_machine_learning_in_r/machine-learning-models.html
# https://dpmartin42.github.io/posts/r/imbalanced-classes-part-1
# https://www.r-bloggers.com/2017/01/handling-class-imbalance-with-r-and-caret-caveats-when-using-the-auc-2/

# Connectivity
# https://cran.r-project.org/web/packages/riverconn/vignettes/Tutorial.html#generalized-riverscape-connectivity-index
# https://cran.r-project.org/web/packages/riverconn/riverconn.pdf
# https://damianobaldan.github.io/riverconn_tutorial/#reach-scale-indices
# https://igraph.org/r/doc/graph_from_data_frame.html
# https://igraph.org/r/doc/set_vertex_attr.html

# --------------------------------------------------------------------- #

# --------------------------------------------------------------------- #
# 00. Set up workspace
# ---------------------------------------------------------------------

# Summary: Set up workspace, load data, and load relevant packages.

# Clean workspace
rm(list = ls())

# Load Packages
if(!require("rstudioapi")){
  install.packages("rstudioapi", dependencies = TRUE); library(rstudioapi)}
if(!require("dplyr")){
  install.packages("dplyr", dependencies = TRUE); library(dplyr)}
if(!require("tidyr")){
  install.packages("tidyr", dependencies = TRUE); library(tidyr)}
if(!require("lubridate")){
  install.packages("lubridate", dependencies = TRUE); library(lubridate)}
if(!require("data.table")){
  install.packages("data.table", dependencies = TRUE); library(data.table)}
if(!require("rio")){
  install.packages("rio", dependencies = TRUE); library(rio)}
if(!require("FSA")){
  install.packages("FSA", dependencies = TRUE); library(FSA)}
# Plotting
if(!require("ggplot2")){
  install.packages("ggplot2", dependencies = TRUE); library(ggplot2)}
if(!require("ggnetwork")){
  install.packages("ggnetwork", dependencies = TRUE); library(ggnetwork)}
if(!require("cowplot")){
  install.packages("cowplot", dependencies = TRUE); library(cowplot)}
if(!require("viridis")){
  install.packages("viridis", dependencies = TRUE); library(viridis)}
if(!require("forcats")){
  install.packages("forcats", dependencies = TRUE); library(forcats)}
if(!require("extrafont")){
  install.packages("extrafont", dependencies = TRUE); library(extrafont)}
if(!require("scales")){
  install.packages("scales", dependencies = TRUE); library(scales)}
# Classification modeling
if(!require("caret")){
  install.packages("caret", dependencies = TRUE); library(caret)}
if(!require("gbm")){
  install.packages("gbm", dependencies = TRUE); library(gbm)}
if(!require("pROC")){
  install.packages("pROC", dependencies = TRUE); library(pROC)}
# Connectivity
if(!require("igraph")){
  install.packages("igraph", dependencies = TRUE); library(igraph)}
if(!require("riverconn")){
  install.packages("riverconn", dependencies = TRUE); library(riverconn)}

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 00. Initialize functions
# ---------------------------------------------------------------------

# Summary: Initialize internal functions

# FUNCTION: network_check() - Verify network structure such at all vertices are unique
#                             and all confluences are binary.

network_check <- function(inData, From_field, To_field){
  # Check that From_Nodes are all unique
  check.reaches    <- as.data.frame(inData %>% group_by({{From_field}}) %>% filter(n() > 1))
  
  # Check that nodes occur no more than twice in To_Node (binary junctions)
  check.confluences <- as.data.frame(inData %>% group_by({{To_field}}) %>% filter(n() > 2))
  
  # Check that To_Nodes all come from From_Nodes
  # NOTE: For network generation that produces a non-existent terminal node, there
  #       should be one row in data.outlet. For a network without a non-existent
  #       terminal node, data.outlet should be blank.
  check.outlet <- data.stream %>% filter(!{{To_field}} %in% {{From_field}})
  
  # Print error outputs
  # ifelse(nrow(check.reaches)     > 0, print('NETWORK CHECK: Duplicate vertices (From nodes)'),
  # ifelse(nrow(check.confluences) > 0, print('NETWORK CHECK: Nonbinary confluences'),
  # ifelse(nrow(check.outlet)      > 1, print('NETWORK CHECK: Nonbinary terminal reaches'),
  # ifelse(nrow(check.outlet)     == 1, print('NETWORK CHECK: To_node not present in From_node'),
  #                                     print('Network check complete')))))
  if(nrow(check.reaches)     > 0) {print('NETWORK CHECK: Duplicate vertices (From nodes) detected')}
  if(nrow(check.confluences) > 0) {print('NETWORK CHECK: Nonbinary confluences detected')}
  if(nrow(check.outlet)      > 1) {print('NETWORK CHECK: Nonbinary terminal reaches detected')}
  if(nrow(check.outlet)     == 1) {print('NETWORK CHECK: To_node not present in From_node detected')}
  print('Network check complete')
  
  return(list(Duplicate_Reaches     = check.reaches,
              Nonbinary_Confluences = check.confluences,
              Terminal_Reaches      = check.outlet))
}


# FUNCTION: generate_attributed_igraph() - Generate an igraph object with edge and vertex attributes,
#                                         add passability fields, and set graph directionality.

generate_attributed_igraph <- function(inData,
                                       From_field,
                                       To_field,
                                       EdgeType_field,
                                       Edge_attributes,
                                       Node_attributes,
                                       Outlet_node,
                                       graphFile){
  # Select edge data and remove nodes not associated with a reach
  data.edges <- inData %>% select({{From_field}},
                                  {{To_field}},
                                  {{Edge_attributes}}) %>%
    filter({{To_field}} %in% {{From_field}})
  
  # Set edge type field
  data.edges <- data.edges %>% mutate(type = ifelse(get({{EdgeType_field}}) == "",
                                                    'Confluence',
                                                    get({{EdgeType_field}})))
  
  # Generate igraph object with edge attributes
  data.graph <- graph_from_data_frame(data.edges)
  
  # Select node data
  data.nodes <- inData %>% select({{From_field}},
                                  {{Node_attributes}})
  
  # Attribute nodes
  for(col in colnames(data.nodes)){
    data.graph <- set_vertex_attr(data.graph,
                                  name = col,
                                  index = V(data.graph),
                                  value = sapply(V(data.graph)$name, function(x){
                                    unlist(data.nodes %>%
                                             filter(From_Node == x) %>%
                                             .[col])
                                  }))
  }
  
  # Assign network directionality based on outlet reach
  data.graph <- set_graph_directionality(data.graph,
                                         field_name = 'name',
                                         outlet_name = as.character(Outlet_node))
  
  # Initialize passability fields
  field.pass <- c('pass_u', 'pass_d')
  for(i in 1:length(field.pass)){
    data.graph <- set_edge_attr(data.graph,
                                field.pass[i],
                                value = 1.0)
  }
  
  # Identify outlet edge for plotting
  index.outlet <- which(V(data.graph)$name == Outlet_node)
  
  # Set plotting dimensions
  size.plot <- data.frame(node  = NA,
                          edge  = NA,
                          arrow = NA,
                          text  = NA)
  ifelse(length(V(data.graph)) <= 100,  size.plot[1,] <- c(0.1, 1, 10, 5),
         ifelse(length(V(data.graph)) <= 1000, size.plot[1,] <- c(0.05, 0.5, 5, 2.5),
                size.plot[1,] <- c(0.01, 0.1, 1, 0.5)))
  
  # Plot to confirm
  gg0 <- ggnetwork(data.graph,
                   layout =  layout_as_tree(data.graph %>% as.undirected, root = index.outlet),
                   scale = FALSE)
  plot <-
    ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_nodes(alpha = 0.3,
               size = size.plot$node) +
    geom_edges(alpha = 0.5,
               arrow = arrow(length = unit(size.plot$arrow, "pt"), type = "closed"),
               linewidth = size.plot$edge,
               aes(color = type)) +
    scale_color_viridis(discrete = TRUE)+
    geom_nodetext(aes(label = name), fontface = "bold",
                  size = size.plot$text) +
    theme_blank()
  ggsave(graphFile, plot = plot,
         width = 40, height = 40, units = 'cm',
         dpi = 600)
  
  return(data.graph)
  
}


# FUNCTION: get_frequency_summary() - Function that generates a frequency tables from supplied groupings

get_frequency_summary <- function(data, fields.Group){
  data.out <- data %>% group_by_at({{fields.Group}}) %>%
    summarise(n = n()) %>%
    mutate(freq = round(n / sum(n), digits = 2)) %>%
    as.data.frame
  return(data.out)
}


# FUNCTION: generate_random_pass() - Function that randomly assigns passability ratings to barriers
#                                    given an expected barrier type and observed probability

generate_random_pass <- function(data, # Dataframe of barriers to predict
                                 freq.tbl, # Table of frequencies used for random sampling
                                 n, # Number of random samples to perform
                                 seed.start, # Seed to start random samples
                                 field.Type, # Field that maps the barrier type frequency
                                 field.Pass, # Field that maps the barrier passability
                                 field.Predict,
                                 field.ID){ # Field used to separate which records should be predicted
  
  # Initialize output list
  out.list <- list()
  
  # Generate seeds for reproducible random assignments
  set.seed(seed.start)
  seed.Random <- sample(x = 1000000, size = n, replace = FALSE)
  
  # Loop across all random assignment seeds
  for(i in 1:length(seed.Random)){
    # Initialize output data
    data.out <- data %>% mutate(Pass_R = !!sym({{field.Pass}}))
    
    # Select data for prediction
    data.predict <- data.out %>% filter(!!sym({{field.Predict}}) == "" | is.na(!!sym({{field.Predict}})))
    
    # Extract type field values to predict from 
    types.Pred <- data.predict %>% distinct(!!sym({{field.Type}})) %>% pull(!!sym({{field.Type}}))
    
    # Assign random passabilities by looping across each type
    for(j in 1:length(types.Pred)){
      # Extract data by type
      data.type <- data.predict %>% filter(!!sym({{field.Type}}) == types.Pred[j])
      
      # Extract the passability ratings for a given type
      type <- freq.tbl %>% filter(!!sym({{field.Type}}) == types.Pred[j]) %>% pull(!!sym({{field.Pass}}))
      
      # Extract passability rating probabilities for given type
      prob.type <- freq.tbl %>% filter(!!sym({{field.Type}}) == types.Pred[j]) %>% pull(freq)
      
      # Randomly assign values to dataset
      set.seed(seed.Random[i])
      data.type <- data.type %>% mutate(Pass_R = sample(x = type,
                                                        size = nrow(data.type),
                                                        prob = prob.type,
                                                        replace = TRUE)) %>%
        select(!!sym({{field.ID}}), Pass_R)
      
      # Join data to output table
      data.out <- data.out %>% left_join(data.type, by = as.character({{field.ID}})) %>%
        mutate(Pass_R = coalesce(Pass_R.y, Pass_R.x)) %>%
        select(-Pass_R.y, -Pass_R.x)
    }
    
    # Name output as seed
    name.Seed <- paste0('Seed_', as.character(seed.Random[i]))
    
    # Assign barrier randomization to output list
    out.list[[name.Seed]] <- data.out
  }
  
  # Return output data
  return(out.list)
}


# FUNCTION: generate_binary_model() - Function that generates a binary gradient boosting classification
#                                     model for input barrier presence/absence from predictor variables.

generate_binary_model <- function(data.sample,
                                  response,
                                  predictors,
                                  seed.split,
                                  perc.train,
                                  n.fold,
                                  type.sample,
                                  type.model,
                                  type.metric){
  # Initialize output list
  out.list <- list()
  
  # Select data for model
  data.model <- data.sample %>% select(all_of({{response}}), all_of({{predictors}})) %>%
    mutate_if(is.character, as.factor)
  
  # Split train/test datasets
  set.seed(seed.split)
  train.Index <- createDataPartition(data.model[,{{response}}],
                                     p = perc.train,
                                     list = FALSE, times = 1)
  
  # Generate train/test split
  data.train <- data.model[train.Index,]
  data.test <- data.model[-train.Index,]
  
  # Tune model parameters
  ctrl <- trainControl(method = 'cv', number = n.fold, returnResamp = 'none',
                       summaryFunction = twoClassSummary, classProbs = TRUE)
  
  # Downsample for model training
  ctrl$sampling <- type.sample
  
  # Declare formula
  model.formula <- reformulate('.', response = response)
  
  # Train model
  model.out <- train(model.formula,
                     data = data.train,
                     method = type.model,
                     trControl = ctrl,
                     metric = type.metric, 
                     verbose = FALSE)
  
  # Model evaluation - accuracy and kappa
  predictions <- predict(object = model.out, data.test[,{{predictors}}], type = 'raw')
  predictions.cm <- confusionMatrix(data.test[,{{response}}], as.factor(predictions))
  # print(predictions.cm$overall[1])
  
  # Model evaluation - AUC
  predictions.prob <- predict(object = model.out, data.test[,{{predictors}}], type = 'prob')
  predictions.auc <- roc(ifelse(data.test[,{{response}}] == 'Present', 1, 0), predictions.prob[[2]])
  # print(predictions.auc$auc)
  
  # Generate model evaluation dataframe
  tbl.eval <- data.frame(Model_seed = seed.split,
                         Accuracy = predictions.cm$overall[1],
                         Kappa = predictions.cm$overall[2],
                         AUC = predictions.auc$auc[1], 
                         row.names = NULL)
  
  # Generate names for list output
  name.model <- paste0('Model_Seed_', seed.split)
  name.eval  <- paste0('Eval_Seed_', seed.split)
  
  # Assign outputs to list
  out.list[[name.model]] <- model.out
  out.list[[name.eval]] <- tbl.eval
  
  # Return output model
  return(out.list)
}


# FUNCTION: generate_multiclass_model() - Function that generates a binary gradient boosting classification
#                                         model for input barrier passability from predictor variables. 

generate_multiclass_model <- function(data.sample,
                                      response,
                                      predictors,
                                      seed.split,
                                      perc.train,
                                      n.fold,
                                      type.sample,
                                      type.model,
                                      type.metric){
  # Initialize output list
  out.list <- list()
  
  # Select data for model
  data.model <- data.sample %>% select(all_of({{response}}), all_of({{predictors}})) %>%
    mutate_if(is.character, as.factor)
  
  # Split train/test datasets
  set.seed(seed.split)
  train.Index <- createDataPartition(data.model[,{{response}}],
                                     p = perc.train,
                                     list = FALSE, times = 1)
  
  # Generate train/test split
  data.train <- data.model[train.Index,]
  data.test <- data.model[-train.Index,]
  
  # Tune model parameters
  ctrl <- trainControl(method = 'cv', number = n.fold, returnResamp = 'none',
                       summaryFunction = multiClassSummary, classProbs = TRUE)
  
  # Downsample for model training
  ctrl$sampling <- type.sample
  
  # Declare formula
  model.formula <- reformulate('.', response = response)
  
  # Train model
  model.out <- train(model.formula,
                     data = data.train,
                     method = type.model,
                     trControl = ctrl,
                     metric = type.metric,
                     verbose = FALSE)
  
  # Model evaluation - accuracy and kappa
  predictions <- predict(object = model.out, data.test[,{{predictors}}], type = 'raw')
  predictions.cm <- confusionMatrix(data.test[,{{response}}], as.factor(predictions))
  # print(predictions.cm$overall[1])
  
  # Model evaluation - AUC
  predictions.prob <- predict(object = model.out, data.test[,{{predictors}}], type = 'prob')
  predictions.auc <- multiclass.roc(data.test[,{{response}}], predictions.prob)
  # print(predictions.auc$auc)
  
  # Generate model evaluation dataframe
  tbl.eval <- data.frame(Model_seed = seed.split,
                         Accuracy = predictions.cm$overall[1],
                         Kappa = predictions.cm$overall[2],
                         AUC = predictions.auc$auc[1], 
                         row.names = NULL)
  
  # Generate names for list output
  name.model <- paste0('Model_Seed_', seed.split)
  name.eval  <- paste0('Eval_Seed_', seed.split)
  
  # Assign outputs to list
  out.list[[name.model]] <- model.out
  out.list[[name.eval]] <- tbl.eval
  
  # Return output model
  return(out.list)
}


# FUNCTION: generate_modeled_barriers() - Function that fits binary or multiclass gradient boosting
#                                         classification models to predict barrier presence/absence
#                                         or passability from predictors, then predicts barrier passability
#                                         for unsampled barriers in a prediction dataset.

generate_modeled_barriers <- function(data.sample,
                                      data.predict,
                                      type.predicted,
                                      type.observed,
                                      pass.observed,
                                      field.id,
                                      response,
                                      predictors,
                                      seed.start,
                                      perc.train,
                                      n.fold,
                                      n.models,
                                      type.sample,
                                      type.model,
                                      type.metric,
                                      str.model){
  # Initialize output list
  list.out <- list()
  
  # Generate seeds
  set.seed(seed.start)
  seed.Model <- sample(x = 1000000, size = n.models, replace = FALSE)
  # print(paste0(c('Seed.model: ', seed.Model), collapse = " "))
  
  # Separate loops based on binary/multiclass model
  if (str.model == 'binary') {
    print('Binary model')
    for(i in 1:length(seed.Model)){
      # Generate binary model
      model.out <- generate_binary_model(data.sample = data.sample,
                                         response = response,
                                         predictors = predictors,
                                         seed.split = seed.Model[i],
                                         perc.train = perc.train,
                                         n.fold = n.fold,
                                         type.sample = type.sample,
                                         type.model = type.model,
                                         type.metric = type.metric)
      
      # Generate name assignments model and evaluation for output
      name.out <- paste0('binary_seed_', seed.Model[i])
      
      # Attach model and evaluation output to list
      list.out[['Model']][[name.out]] <- model.out[[1]]
      list.out[['Eval']][[name.out]] <- model.out[[2]]
      
      # Declare output data
      barriers.out <- data.predict %>% select(all_of({{field.id}}), 
                                              all_of({{type.predicted}}), 
                                              all_of({{type.observed}}),
                                              all_of({{pass.observed}}),
                                              all_of({{predictors}}))
      
      # Initialize predicted passability fields
      barriers.out <- barriers.out %>% mutate(Model_val = NA,
                                              Pass_M = NA)
      
      # Assing known passability ratings to model prediction
      barriers.out <- barriers.out %>% mutate(Pass_M = !!sym({{pass.observed}}))
      
      # Declare data for model prediction
      barriers.predict <- barriers.out %>% filter(is.na(Pass_M))
      
      # Predict barriers with model
      barriers.predict$Model_val <- predict(object = model.out[[1]],
                                            barriers.predict[,{{predictors}}],
                                            type = 'raw')
      
      # Assign passability based on barrier presence
      barriers.predict <- barriers.predict %>% mutate(Pass_M = ifelse(Model_val == 'Absent', '1', '0.5'))
      
      # Reduce join to only related data
      barriers.join <- barriers.predict %>% select(all_of({{field.id}}), Model_val, Pass_M)
      # print(head(barriers.join))
      
      # Join predictions to output barriers
      barriers.out <- barriers.out %>% left_join(barriers.join, by = field.id) %>%
        mutate(Model_val = if_else(is.na(Model_val.x), Model_val.y, Model_val.x),
               Pass_M = ifelse(is.na(Pass_M.x), Pass_M.y, Pass_M.x)) %>%
        select(-Model_val.x, -Model_val.y, -Pass_M.x, -Pass_M.y)
      
      # Attach predicted barriers to list
      list.out[['Model_Results']][[name.out]] <- barriers.predict
      list.out[['Barriers']][[name.out]] <- barriers.out
      
      # Print completion log
      print(paste0('COMPLETE: ', as.character(i), ' / ', as.character(length(seed.Model))))
      
    }
    
  } else if (str.model == 'multiclass'){
    print('Multiclass model')
    for(i in 1:length(seed.Model)){
      # Generate binary model
      model.out <- generate_multiclass_model(data.sample = data.sample,
                                             response = response,
                                             predictors = predictors,
                                             seed.split = seed.Model[i],
                                             perc.train = perc.train,
                                             n.fold = n.fold,
                                             type.sample = type.sample,
                                             type.model = type.model,
                                             type.metric = type.metric)
      
      # Generate name assignments model and evaluation for output
      name.out <- paste0('binary_seed_', seed.Model[i])
      
      # Attach model and evaluation output to list
      list.out[['Model']][[name.out]] <- model.out[[1]]
      list.out[['Eval']][[name.out]] <- model.out[[2]]
      
      # Declare output data
      barriers.out <- data.predict %>% select(all_of({{field.id}}), 
                                              all_of({{type.predicted}}), 
                                              all_of({{type.observed}}),
                                              all_of({{pass.observed}}),
                                              all_of({{predictors}}))
      
      # Initialize predicted passability fields
      barriers.out <- barriers.out %>% mutate(Model_val = NA,
                                              Pass_M = NA)
      
      # Assing known passability ratings to model prediction
      barriers.out <- barriers.out %>% mutate(Pass_M = !!sym({{pass.observed}}))
      
      # Declare data for model prediction
      barriers.predict <- barriers.out %>% filter(is.na(Pass_M))
      
      # Predict barriers with model
      barriers.predict$Model_val <- predict(object = model.out[[1]],
                                            barriers.predict[,{{predictors}}],
                                            type = 'raw')
      
      # Assign passability based on barrier presence
      barriers.predict <- barriers.predict %>% mutate(Pass_M = substring(Model_val, 3))
      
      # Reduce join to only related data
      barriers.join <- barriers.predict %>% select(all_of({{field.id}}), Model_val, Pass_M)
      # print(head(barriers.join))
      
      # Join predictions to output barriers
      barriers.out <- barriers.out %>% left_join(barriers.join, by = field.id) %>%
        mutate(Model_val = if_else(is.na(Model_val.x), Model_val.y, Model_val.x),
               Pass_M = ifelse(is.na(Pass_M.x), Pass_M.y, Pass_M.x)) %>%
        select(-Model_val.x, -Model_val.y, -Pass_M.x, -Pass_M.y)
      
      # Attach predicted barriers to list
      list.out[['Model_Results']][[name.out]] <- barriers.predict
      list.out[['Barriers']][[name.out]] <- barriers.out
      
      # Print completion log
      print(paste0('COMPLETE: ', as.character(i), ' / ', as.character(length(seed.Model))))
    }
    
  } else {
    print('Unsupported model type')
  }
  
  # Return output
  return(list.out)
}


# FUNCTION: join_edge_attributes() - Function that joins edge attributes based on a common identifier field.
#                                    Function currently hard-coded for 'UID', 'Type_R', and 'Pass_R' fields.

join_edge_attributes <- function(inGraph, inData, field.Pass){
  # Initialize output
  outGraph <- inGraph
  
  # Join pass field
  outGraph <- set_edge_attr(outGraph,
                            name = as.name(field.Pass),
                            index = E(outGraph),
                            value = as.numeric(sapply(E(outGraph)$UID, function(x){
                              unlist(inData %>%
                                       filter(UID == x) %>%
                                       .[[as.name(field.Pass)]])
                            })))
  
  # Convert pass field to numeric
  outGraph <- set_edge_attr(outGraph,
                            name = as.name(field.Pass),
                            index = E(outGraph),
                            value = ifelse(E(outGraph)$type == 'Confluence',
                                           1.0, E(outGraph)$Pass_R))
  
  # Return output
  return(outGraph)
}

# FUNCTION: calculate_random_dci() - Function calculates symmetric and aysmmetric DCI from a provided igraph object,
#                                    dataframe of barrier passability, and fields indicating the
#                                    passability, and weight fields of the input tables.

calculate_dci <- function(inGraph, inBarriers, scenario, field.Pass, field.Weight, field.ID){
  # Initialize output table
  data.out <- data.frame(Scenario = character(0),
                         ID = character(0),
                         DCI_symm = numeric(0),
                         DCI_asym = numeric(0))
  
  # Join barrier data to dataframe
  graph.out <- join_edge_attributes(inGraph = inGraph,
                                    inData = inBarriers,
                                    field.Pass = field.Pass)
  
  # Attribute pass_u and pass_d with joined passability values
  field.pass <- c('pass_u', 'pass_d')
  for(i in 1:length(field.pass)){
    graph.out <- set_edge_attr(graph.out,
                               field.pass[i],
                               value = E(graph.out)$Pass_R)
  }
  
  # Calculate DCI.symmetric
  dci.symm <- index_calculation(graph = graph.out,
                                weight = field.Weight,
                                B_ij_flag = FALSE,
                                index_type = 'full',
                                dir_fragmentation_type = 'symmetric')
  
  # Calculate DCI.asymmetric
  dci.asym <- index_calculation(graph = graph.out,
                                weight = field.Weight,
                                B_ij_flag = FALSE,
                                index_type = 'full',
                                dir_fragmentation_type = 'asymmetric')
  
  # Attribute output
  data.out <- data.out %>% add_row(Scenario = scenario,
                                   ID = field.ID,
                                   DCI_symm = dci.symm$index,
                                   DCI_asym = dci.asym$index)
  
  # Return output
  # return(graph.out)
  return(data.out)
}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.2. Barrier data and detection
# ---------------------------------------------------------------------

# Set inputs
raw.data            <- 'Stream_BR.csv'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv(raw.data,
                        header = TRUE)

# Sampled barriers
data.barriers <- data.stream %>% filter(UID != "")

# Summarize barrier survey by sournce
source.sum <- data.barriers %>% group_by(SOURCE, Barrier_Expected) %>%
  summarise(n = n()) %>% as.data.frame

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.5.a. Barrier modeling - Declare data
# ---------------------------------------------------------------------

# Summary: This section declares data for the presence/absence, presence/absence,
#          and passability barrier scenario models.

# Set inputs
raw.data            <- 'Stream_BR.csv'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv(raw.data,
                        header = TRUE)

# Identify barriers
data.barriers <- data.stream %>% filter(UID != "")

# Identify sampled barriers
data.sample <- data.stream %>% filter(Sampled == TRUE & Pass != '')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.5.b. Barrier modeling - Uniform method method
# ---------------------------------------------------------------------

# Summary: This section generates a barrier dataset with uniform passability ratings 

# Initialize output list
barriers.uniform <- list()

# Define range of uniform passability values
pass.range <- seq(0, 0.99, by = 0.01)

# Uniformly assign passability ratings across range
for(i in 1:length(pass.range)){
  # Name scenario
  seed.id <- paste0('UniPass_', as.character(pass.range[i]))

  # Select barriers
  barriers.out <- data.barriers %>% select(UID, Barrier_Expected, Pass)
  
  # Add uniform passability rating and rating ID
  barriers.out <- barriers.out %>% mutate(Pass_R = ifelse((is.na(Pass) | Pass == ""),
                                                          pass.range[i],
                                                          Pass),
                                          Uniform_ID = seed.id)
  
  # Add to output list
  barriers.uniform[[seed.id]] <- barriers.out}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.5.c. Barrier modeling - Presence/absence method
# ---------------------------------------------------------------------

# Summary: This section generates the presence/absence barrier models.

# Add barrier presence/absence field
data.sample.model <- data.sample %>% mutate(Pass_chr = paste0('P_', as.character(Pass)),
                                            Barrier_Present = ifelse(Pass_chr == 'P_1', 'Absent', 'Present'))

# Log10(n+1) transform continuous predictors
data.sample.model <- data.sample.model %>% mutate_at(.vars = c('Elev_m',
                                                                'Slope_site_perc',
                                                                'Slope_reach_perc',
                                                                'Slope_segment_perc',
                                                                'Qmad_cfs'),
                                                     .funs = function(x) log10(x + 1))

# Set model inputs
model.data.sample    <- data.sample.model
model.response       <- 'Barrier_Present'
model.predictors     <- c('tnmfrc', 'Elev_m', 'Slope_site_perc',
                          'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')

model.type.sample    <- 'down'
model.type.model     <- 'gbm'
model.type.metric    <- 'ROC'
model.perc.train     <- 0.7
model.n.fold         <- 5

model.type.observed  <- 'Barrier_Observed'
model.type.predicted <- 'Barrier_Expected'
model.pass.observed  <- 'Pass'
model.field.id       <- 'UID'

model.data.predict   <- data.barriers %>% select(all_of({{model.field.id}}),
                                                 all_of({{model.type.observed}}),
                                                 all_of({{model.type.predicted}}),
                                                 all_of({{model.pass.observed}}),
                                                 all_of({{model.predictors}}))
model.seed.start     <- 223445
model.n.models       <- 100
model.str.model      <-'binary'


# Generate binary barrier model
barriers.presabs <- generate_modeled_barriers(data.sample = data.sample.model,
                                              response = model.response,
                                              predictors = model.predictors,
                                              type.sample = model.type.sample,
                                              type.model = model.type.model,
                                              type.metric = model.type.metric,
                                              perc.train = model.perc.train,
                                              n.fold = model.n.fold,
                                              
                                              type.observed = model.type.observed,
                                              type.predicted = model.type.predicted,
                                              pass.observed = model.pass.observed,
                                              field.id = model.field.id,
                                              
                                              data.predict = model.data.predict,
                                              seed.start = model.seed.start,
                                              n.models = model.n.models,
                                              str.model = model.str.model)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.5.d. Barrier modeling - Rating category method
# ---------------------------------------------------------------------

# Summary: This section generates the barrier passability rating models.

# Add barrier presence/absence field
data.sample.model <- data.sample %>% mutate(Pass_chr = paste0('P_', as.character(Pass)),
                                            Barrier_Present = ifelse(Pass_chr == 'P_1', 'Absent', 'Present'))

# Log10(n+1) transform continuous predictors
data.sample.model <- data.sample.model %>% mutate_at(.vars = c('Elev_m',
                                                               'Slope_site_perc',
                                                               'Slope_reach_perc',
                                                               'Slope_segment_perc',
                                                               'Qmad_cfs'),
                                                     .funs = function(x) log10(x + 1))

# Set model inputs
model.data.sample    <- data.sample.model
model.response       <- 'Pass_chr'
model.predictors     <- c('tnmfrc', 'Elev_m', 'Slope_site_perc',
                          'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')

model.type.sample    <- 'up'
model.type.model     <- 'gbm'
model.type.metric    <- 'ROC'
model.perc.train     <- 0.7
model.n.fold         <- 5

model.type.observed  <- 'Barrier_Observed'
model.type.predicted <- 'Barrier_Expected'
model.pass.observed  <- 'Pass'
model.field.id       <- 'UID'

model.data.predict   <- data.barriers %>% select(all_of({{model.field.id}}),
                                                 all_of({{model.type.observed}}),
                                                 all_of({{model.type.predicted}}),
                                                 all_of({{model.pass.observed}}),
                                                 all_of({{model.predictors}}))
model.seed.start     <- 223445
model.n.models       <- 100
model.str.model      <-'multiclass'

# Generate binary barrier model
barriers.pass <- generate_modeled_barriers(data.sample = data.sample.model,
                                           response = model.response,
                                           predictors = model.predictors,
                                           type.sample = model.type.sample,
                                           type.model = model.type.model,
                                           type.metric = model.type.metric,
                                           perc.train = model.perc.train,
                                           n.fold = model.n.fold,
                                              
                                           type.observed = model.type.observed,
                                           type.predicted = model.type.predicted,
                                           pass.observed = model.pass.observed,
                                           field.id = model.field.id,
                                              
                                           data.predict = model.data.predict,
                                           seed.start = model.seed.start,
                                           n.models = model.n.models,
                                           str.model = model.str.model)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.5.e. Barrier modeling - Random method
# ---------------------------------------------------------------------

# Summary: This section generates the random type and passability assignment models.

# Generate frequency table for expected barrier type
freq.Type <- get_frequency_summary(data = data.sample, fields.Group = c('Barrier_Expected', 'Pass'))


# Generate random barrier type and passability rating assignments from observed distribution
barriers.random <- generate_random_pass(data = data.barriers,
                                        freq.tbl = freq.Type,
                                        seed.start = 223445,
                                        n = 1000,
                                        field.Type = 'Barrier_Expected',
                                        field.Pass = 'Pass',
                                        field.Predict = 'Pass',
                                        field.ID = 'UID')

# Evaluate output type and pass frequencies
# data.check <- barriers.random$Seed_124413
# 
# data.check %>% group_by(Barrier_Expected, Pass_R) %>%
#   summarise(n = n()) %>%
#   mutate(freq = round(n / sum(n), digits = 2)) %>%
#   as.data.frame
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.a Connectivity - Load network and generate igraph object
# ---------------------------------------------------------------------

# Summary: Load stream network, convert to igraph, and calculate no-barrier
#          DCI for validation.

# Set inputs
raw.data            <- 'Stream_BR.csv'

graph.fields        <- c('From_Node', 'To_Node', 'Length_km', 'UID', 'Barrier_Expected') # Fields to build graph with
attributes.edgeType <- 'Barrier_Expected' # String corresponding to barrier type. 
attributes.edge     <- c('UID', 'Barrier_Expected') # Strings corresponding to edge (barrier) attributes.
attributes.node     <- c('From_Node', 'Length_km') # String corresponding to node (stream segment) attributes. A field to weight connectivity indices is required, commonly assessed as length or HSI.
field.weight        <- 'Length_km' # String corresponding to node attribute used to weight connectivity indices.
file.graph          <- 'Graph_BR.png' # String indicating the name of the output file used to check igraph structure.

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv(raw.data,
                        header = TRUE)

# Filter data for igraph
# NOTE: Input data represents the edge list
data.graph <- data.stream %>% select({{graph.fields}})

# Check network structure
data.check <- network_check(inData = data.graph,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node as the 'From_Node' on the terminal reach
outlet <- data.check$Terminal_Reaches$From_Node

# Generate igraph object
graph.stream <- generate_attributed_igraph(inData          = data.stream,
                                           From_field      = From_Node,
                                           To_field        = To_Node,
                                           EdgeType_field  = {{attributes.edgeType}},
                                           Edge_attributes = {{attributes.edge}},
                                           Node_attributes = {{attributes.node}},
                                           Outlet_node     = outlet,
                                           graphFile       = file.graph)

# Validate igraph object - All pass_u and pass_d are set to 1 (fully passable),
#                          so DCI should equal 0 (i.e. no fragementation)
index_calculation(graph = graph.stream,
                  weight = field.weight,
                  B_ij_flag = FALSE,
                  dir_fragmentation_type = 'symmetric')

# View(as_data_frame(graph.stream, what = 'edges'))
# View(as_data_frame(graph.stream, what = 'vertices'))

# # Initialize results output
# results.connectivity <- list('Date_Generated' = Sys.time())
# 
# # Export list
# export(results.connectivity, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.b Connectivity - Natural barriers
# ---------------------------------------------------------------------

# Summary: Calculate DCI for natural barriers

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Set inputs
outScenario <- 'Natural'

# Initialize output
connectivity.out <- data.frame(Scenario = outScenario,
                               Seed = NA,
                               DCI_symm = NA,
                               DCI_asym = NA)

# Set all natural barriers to impassable (pass_u = pass_d = 0.0)
field.pass <- c('pass_u', 'pass_d')
graph.conn <- graph.stream
for(i in 1:length(field.pass)){
  graph.conn <- set_edge_attr(graph.conn,
                              field.pass[i],
                              value = ifelse(E(graph.conn)$type == 'Waterfall',
                                             0.0, 1.0))}

# View(as_data_frame(graph.conn, what = 'edges'))
# View(as_data_frame(graph.conn, what = 'vertices'))

# Calculate connectivity indices
conn.dci_symm <- index_calculation(graph = graph.conn,
                                   weight = field.weight,
                                   index_type = 'full',
                                   B_ij_flag = FALSE,
                                   dir_fragmentation_type = 'symmetric')

conn.dci_asym <- index_calculation(graph = graph.conn,
                                   weight = field.weight,
                                   B_ij_flag = FALSE,
                                   index_type = 'full',
                                   dir_fragmentation_type = 'asymmetric')

# Add DCI to output table
connectivity.out <- connectivity.out %>% mutate(DCI_symm = conn.dci_symm$index,
                                                DCI_asym = conn.dci_asym$index)

# Assign output table to connectivity output
results.connectivity[[outScenario]] <- connectivity.out

# Save output
export(results.connectivity, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.c Connectivity - Dam barriers
# ---------------------------------------------------------------------

# Summary: Calculate DCI for dam barriers set to impassable (pass_u = pass_d = 0.0)

# Set inputs
outScenario <- 'Dam'

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Initialize output
connectivity.out <- data.frame(Scenario = outScenario,
                               Seed = NA,
                               DCI_symm = NA,
                               DCI_asym = NA)

# Set all natural barriers to impassable (pass_u = pass_d = 0.0)
field.pass <- c('pass_u', 'pass_d')
graph.conn <- graph.stream
for(i in 1:length(field.pass)){
  graph.conn <- set_edge_attr(graph.conn,
                              field.pass[i],
                              value = ifelse((E(graph.conn)$type == 'Waterfall') | (E(graph.conn)$type == 'Dam') ,
                                             0.0, 1.0))}

# View(as_data_frame(graph.conn, what = 'edges'))
# View(as_data_frame(graph.conn, what = 'vertices'))

# Calculate connectivity indices
conn.dci_symm <- index_calculation(graph = graph.conn,
                                   weight = field.weight,
                                   index_type = 'full',
                                   B_ij_flag = FALSE,
                                   dir_fragmentation_type = 'symmetric')

conn.dci_asym <- index_calculation(graph = graph.conn,
                                   weight = field.weight,
                                   B_ij_flag = FALSE,
                                   index_type = 'full',
                                   dir_fragmentation_type = 'asymmetric')

# Add DCI to output table
connectivity.out <- connectivity.out %>% mutate(DCI_symm = conn.dci_symm$index,
                                                DCI_asym = conn.dci_asym$index)

# Assign output table to connectivity output
results.connectivity[[outScenario]] <- connectivity.out

# Save output
export(results.connectivity, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.d Connectivity - Uniform method
# ---------------------------------------------------------------------

# Summary: Calculate DCI with uniform barrier passability rating assignments

# Set inputs
connectivity.barriers <- barriers.uniform
connectivity.scenario     <- 'Uniform'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Initialize output
connectivity.out <- data.frame(Scenario = character(0),
                               Seed = character(0),
                               DCI_symm = numeric(0),
                               DCI_asym = numeric(0))

# Calculate DCI and add results to output across barrier randomizations
for(i in 1:length(connectivity.barriers)){
  
  # Name scenario for output
  seed.id <- names(connectivity.barriers[i])
  
  # Calculate DCI
  dci.out <- calculate_dci(inGraph = graph.stream,
                           inBarriers = connectivity.barriers[[i]],
                           scenario = connectivity.scenario,
                           field.Pass = connectivity.field.pass,
                           field.Weight = connectivity.field.weight,
                           field.ID = seed.id)
  
  # Rename ID field as Seed
  dci.out <- dci.out %>% rename(Seed = ID) %>% as.data.frame
  
  # Add DCI calculate to output
  connectivity.out <- rbind(connectivity.out, dci.out)
  
  # Print completion  log
  print(paste0('COMPLETE: ', as.character(i), ' / ', as.character(length(connectivity.barriers))))}

# Assign output table to connectivity output
results.connectivity[[connectivity.scenario]] <- connectivity.out

# Save output
export(results.connectivity, file = 'Results_Connectivity.xlsx')


# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.e Connectivity - Random method
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Set inputs
connectivity.barriers <- barriers.random
connectivity.scenario     <- 'Random'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Initialize output
connectivity.out <- data.frame(Scenario = character(0),
                               Seed = character(0),
                               DCI_symm = numeric(0),
                               DCI_asym = numeric(0))

# Calculate DCI and add results to output across barrier randomizations
for(i in 1:length(connectivity.barriers)){
  
  # Name scenario for output
  seed.id <- names(connectivity.barriers[i])
  
  # Calculate DCI
  dci.out <- calculate_dci(inGraph = graph.stream,
                           inBarriers = connectivity.barriers[[i]],
                           scenario = connectivity.scenario,
                           field.Pass = connectivity.field.pass,
                           field.Weight = connectivity.field.weight,
                           field.ID = seed.id)
  
  # Rename ID field as Seed
  dci.out <- dci.out %>% rename(Seed = ID) %>% as.data.frame
  
  # Add DCI calculate to output
  connectivity.out <- rbind(connectivity.out, dci.out)
  
  # Print completion  log
  print(paste0('COMPLETE: ', as.character(i), ' / ', as.character(length(connectivity.barriers))))}

# Assign output table to connectivity output
results.connectivity[[connectivity.scenario]] <- connectivity.out

# Save output
export(results.connectivity, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.f Connectivity - Presence/absence method
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Set inputs
connectivity.barriers <- barriers.presabs$Barriers
connectivity.scenario     <- 'PresAbs'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

i <- 1

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Initialize output
connectivity.out <- data.frame(Scenario = character(0),
                               Seed = character(0),
                               DCI_symm = numeric(0),
                               DCI_asym = numeric(0))

# Calculate DCI and add results to output across barrier randomizations
for(i in 1:length(connectivity.barriers)){
  
  # Generate Pass_R field for barrier passability (required for hard coding in join_edge_attributes)
  connectivity.barriers[[i]][[connectivity.field.pass]] <- connectivity.barriers[[i]][['Pass_M']]
  
  # Convert passability field to numeric
  connectivity.barriers[[i]][[connectivity.field.pass]] <- as.numeric(connectivity.barriers[[i]][[connectivity.field.pass]])
  
  # Name scenario for output
  seed.id <- names(connectivity.barriers[i])
  
  # Calculate DCI
  dci.out <- calculate_dci(inGraph = graph.stream,
                           inBarriers = connectivity.barriers[[i]],
                           scenario = connectivity.scenario,
                           field.Pass = connectivity.field.pass,
                           field.Weight = connectivity.field.weight,
                           field.ID = seed.id)
  
  # Rename ID field as Seed
  dci.out <- dci.out %>% rename(Seed = ID) %>% as.data.frame
  
  # Add DCI calculate to output
  connectivity.out <- rbind(connectivity.out, dci.out)
  
  # Print completion  log
  print(paste0('COMPLETE: ', as.character(i), ' / ', as.character(length(connectivity.barriers))))}

# Assign output table to connectivity output
results.connectivity[[connectivity.scenario]] <- connectivity.out

# Save output
export(results.connectivity, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: 2.6.g Connectivity - Passability rating method
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Set inputs
connectivity.barriers <- barriers.pass$Barriers
connectivity.scenario     <- 'PassRate'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

i <- 1

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Initialize output
connectivity.out <- data.frame(Scenario = character(0),
                               Seed = character(0),
                               DCI_symm = numeric(0),
                               DCI_asym = numeric(0))

# Calculate DCI and add results to output across barrier randomizations
for(i in 1:length(connectivity.barriers)){
  
  # Generate Pass_R field for barrier passability (required for hard coding in join_edge_attributes)
  connectivity.barriers[[i]][[connectivity.field.pass]] <- connectivity.barriers[[i]][['Pass_M']]
  
  # Convert passability field to numeric
  connectivity.barriers[[i]][[connectivity.field.pass]] <- as.numeric(connectivity.barriers[[i]][[connectivity.field.pass]])
  
  # Name scenario for output
  seed.id <- names(connectivity.barriers[i])
  
  # Calculate DCI
  dci.out <- calculate_dci(inGraph = graph.stream,
                           inBarriers = connectivity.barriers[[i]],
                           scenario = connectivity.scenario,
                           field.Pass = connectivity.field.pass,
                           field.Weight = connectivity.field.weight,
                           field.ID = seed.id)
  
  # Rename ID field as Seed
  dci.out <- dci.out %>% rename(Seed = ID) %>% as.data.frame
  
  # Add DCI calculate to output
  connectivity.out <- rbind(connectivity.out, dci.out)
  
  # Print completion  log
  print(paste0('COMPLETE: ', as.character(i), ' / ', as.character(length(connectivity.barriers))))}

# Assign output table to connectivity output
results.connectivity[[connectivity.scenario]] <- connectivity.out

# Save output
export(results.connectivity, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: 3.1. Barrier Surveys
# ---------------------------------------------------------------------

# Set inputs
raw.data  <- 'Stream_BR.csv'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv(raw.data,
                        header = TRUE)

# Sampled barriers
data.sample <- data.stream %>% filter(Sampled == TRUE & Pass != '')


# Summarize barrier survey by type
type.sum <- data.sample %>% group_by(Basin, Barrier_Expected, Barrier_Observed) %>%
                            summarise(n = n()) %>%
                            mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame

# Summarize barrier survey by passability
pass.sum <- data.sample %>% group_by(Basin, Barrier_Observed, Pass) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame

# Summarize barrier survey by all fields
pass.sum <- data.sample %>% group_by(Basin, Barrier_Expected, Barrier_Observed, Pass) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame


# Summarize by barrier observed
type.sum <- data.sample %>% group_by(Barrier_Observed) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame

# Summarize by passability
pass.sum <- data.sample %>% group_by(Barrier_Observed, Pass) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame

# Summarize barrier survey by passability limit
limit.sum <- data.sample %>% group_by(Barrier_Observed, Limit_WDFG) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), digits = 2))

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: 3.2.a. Barrier Modeling - Presence/absence model
# ---------------------------------------------------------------------

# Summary: Summarize performance metrics (auc) and variable importance
#          from gradient boosting models

# Model evaluation - Presence/absence model

# Retrieve model performance tables
eval.presabs <- do.call('rbind', barriers.presabs$Eval)

# Calculate summary statistics for each metric
eval.presabs.sum <- eval.presabs %>% pivot_longer(cols = c('Accuracy', 'Kappa', 'AUC'),
                                                  names_to = 'Metric',
                                                  values_to = 'Value') %>%
                                     group_by(Metric) %>%
                                     summarize(n = n(),
                                               Mean = mean(Value, na.rm = TRUE),
                                               Median = median(Value, na.rm = TRUE),
                                               SD = sd(Value, na.rm = TRUE)) %>% 
                                     mutate(SE = SD / sqrt(n),
                                            ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
                                            ci95.up = Mean + ((qt(0.975, df = n-1)) * SE)) %>%
                                     mutate_if(is.numeric, round, 2) %>%
                                     as.data.frame

# Gather variable importance information for each model
varimp.list <- list()
for(i in 1:length(names(barriers.presabs[[1]]))){
  data.out <- summary.gbm(barriers.presabs[[1]][[i]][[11]])
  varimp.list[[i]] <- data.out
}
eval.presabs.varimp <- as.data.frame(rbindlist(varimp.list, idcol = 'index'))

# Summarize variable importance
eval.presabs.varimp.sum <- eval.presabs.varimp %>% group_by(var) %>%
  summarize(Mean = mean(rel.inf, na.rm = TRUE)) %>%
  arrange(desc(Mean)) %>%
  as.data.frame

# Add ranking for each site
eval.presabs.varimp <- eval.presabs.varimp %>% arrange(index, rel.inf) %>%
                                               group_by(index) %>%
                                               mutate(rank = rank(desc(rel.inf), ties.method = 'max'))

# Summarize variable ranking
eval.presabs.varimp.sum.rank <- eval.presabs.varimp %>% group_by(rank, var) %>%
                                                        summarize(total_count = n(), .groups = 'drop') %>%
                                                        arrange(rank, desc(total_count)) %>%
                                                        as.data.frame

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: 3.2.b. Barrier Modeling - Passability model
# ---------------------------------------------------------------------

# Summary: Summarize performance metrics (auc) and variable importance
#          from gradient boosting models

# Model evaluation - Passability model

# Retrieve model performance tables
eval.pass <- do.call('rbind', barriers.pass$Eval)

# Calculate summary statistics for each metric
eval.pass.sum <- eval.pass %>% pivot_longer(cols = c('Accuracy', 'Kappa', 'AUC'),
                                                      names_to = 'Metric',
                                                      values_to = 'Value') %>%
  group_by(Metric) %>%
  summarize(n = n(),
            Mean = mean(Value, na.rm = TRUE),
            Median = median(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE)) %>% 
  mutate(SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE)) %>%
  mutate_if(is.numeric, round, 2) %>%
  as.data.frame

# Gather variable importance information for each model
varimp.list <- list()
for(i in 1:length(names(barriers.pass[[1]]))){
  data.out <- summary.gbm(barriers.pass[[1]][[i]][[11]])
  varimp.list[[i]] <- data.out
}
eval.pass.varimp <- as.data.frame(rbindlist(varimp.list, idcol = 'index'))

# Add ranking for each site
eval.pass.varimp <- eval.pass.varimp %>% arrange(index, rel.inf) %>%
  group_by(index) %>%
  mutate(rank = rank(desc(rel.inf), ties.method = 'max'))

# Summarize variable ranking
eval.pass.varimp.sum.rank <- eval.pass.varimp %>% group_by(rank, var) %>%
  summarize(total_count = n(), .groups = 'drop') %>%
  arrange(rank, desc(total_count)) %>%
  as.data.frame

# Summarize variable importance
eval.pass.varimp.sum <- eval.pass.varimp %>% group_by(var) %>%
  summarize(Mean = mean(rel.inf, na.rm = TRUE)) %>%
  arrange(desc(Mean)) %>%
  as.data.frame

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: 3.2. Partial dependence plot data
# ---------------------------------------------------------------------

# Summary: Extract information for partial dependence plots for 
#          presence/absence and rating category models.

# Declare input
pred.vars <- c('Elev_m', 'Qmad_cfs', 'Slope_site_perc', 'Slope_reach_perc', 'Slope_segment_perc')
models.presabs <- names(barriers.presabs$Model)
models.pass <- names(barriers.pass$Model)
class.pass <- c('P_0', 'P_0.33', 'P_0.67', 'P_1')

# Initialize output
data.out <- data.frame(Value = numeric(),
                       yhat = numeric(),
                       Model = character(),
                       Variable = character(),
                       Alternative = character(),
                       Class = character())

# Extract partial dependence plot info for each model and continuous predictor
## Presence/absence model
pdp.presabs <- data.out
for(i in 1:length(pred.vars)){
  for(j in 1:length(models.presabs)){
    # Retrieve data
    data.pdp <- pdp::partial(object = barriers.presabs$Model[[j]],
                             pred.var = pred.vars[i],
                             which.class = 'Present',
                             plot = FALSE)
    # Format for output
    pdp.out <- data.pdp %>%
      rename_with(.cols = 1, ~'Value') %>%
      mutate(Model = models.presabs[j],
             Variable = pred.vars[i],
             Alternative = 'Presence/absence',
             Class = 'Present')
    pdp.presabs <- rbind(pdp.presabs, pdp.out)
  }
}

## Rating category model
pdp.pass <- data.out
for(i in 1:length(pred.vars)){
  for(j in 1:length(models.pass)){
    for(k in 1:length(class.pass)){
      # Retrieve data
      data.pdp <- pdp::partial(object = barriers.pass$Model[[j]],
                               pred.var = pred.vars[i],
                               which.class = class.pass[k],
                               plot = FALSE)
      # Format for output
      pdp.out <- data.pdp %>%
        rename_with(.cols = 1, ~'Value') %>%
        mutate(Model = models.presabs[j],
               Variable = pred.vars[i],
               Alternative = 'RatingCategory',
               Class = class.pass[k])
      pdp.pass <- rbind(pdp.pass, pdp.out)
    }
  }
}

# Bind datasets
data.pdp <- rbind(pdp.presabs, pdp.pass)

# Read in data
data.all <- rio::import_list('Results_GBM_Models.xlsx')

# Add to output
data.all[['PartialDependence']] <- data.pdp

# Write output
rio::export(data.all, file = 'Results_GBM_Models.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: 3.2.c. Barrier Modeling - Export results for plotting
# ---------------------------------------------------------------------

# Summary: Summarize performance metrics (auc) and variable importance
#          from gradient boosting models

# Add field to evaluation tables identifying model
eval.presabs <- eval.presabs %>% mutate(Model = 'Presence/Absence')
eval.presabs.varimp <- eval.presabs.varimp %>% mutate(Model = 'Presence/Absence')

eval.pass <- eval.pass %>% mutate(Model = 'Passability')
eval.pass.varimp <- eval.pass.varimp %>% mutate(Model = 'Passability')

# Combine datasets
eval.models <- rbind(eval.presabs, eval.pass)
eval.varimp <- rbind(eval.presabs.varimp, eval.pass.varimp)

# Combine in list
eval.out <- list('Model_Stats' = eval.models,
                 'Variable_Importance' = eval.varimp)

# Write output
export(eval.out, file = 'Results_GBM_Models.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# Results: 3.2.d. - Model accuracy and Variable importance summaries
# ---------------------------------------------------------------------

# Set inputs
data.model <- import_list('Results_GBM_Models.xlsx')

# Separate datasets
data.variables <- data.model$Variable_Importance%>% rename(Method = Model) %>%
  mutate(Method = ifelse(Method == 'Passability', 'Rating category',
                         'Presence/absence'))

# Summarize variable ranking
sum.rank <- data.variables %>% group_by(Method, rank, var) %>%
            summarise(total_count = n(), 
                      max_relinf = max(rel.inf),
                      min_relinf = min(rel.inf),
                      .groups = 'drop') %>%
            arrange(Method, rank, desc(total_count)) %>%
            as.data.frame

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: 3.3.a. Connectivity - Summary and analysis
# ---------------------------------------------------------------------

# Summary: Summarize DCI model info and test for statistical difference

# Set inputs
data.all <- import_list('Results_Connectivity.xlsx')

# Combine distribution datasets
data.dci <- rbind(data.all$Uniform, data.all$PresAbs, data.all$PassRate, data.all$Random)

# Rename model scenarios
data.dci <- data.dci %>% mutate(Method = ifelse(Scenario == 'PresAbs', 'Presence/absence',
                                                ifelse(Scenario == 'PassRate', 'Rating category',
                                                       ifelse(Scenario == 'Random', 'Random', 
                                                              ifelse(Scenario == 'Uniform', 'Uniform', NA)))),
                                DCI_p = DCI_symm * 100)

# Calculate stats
stat.dci <- data.dci %>% group_by(Method) %>%
  summarize(n = n(),
            Mean = mean(DCI_p, na.rm = TRUE),
            Median = median(DCI_p, na.rm = TRUE),
            SD = sd(DCI_p, na.rm = TRUE),
            Min = min(DCI_p),
            Max = max(DCI_p)) %>% 
  mutate(Range = Max - Min,
          SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE),
         Min_Perc = Min / 12.2,
         Max_Perc = Max / 12.2,
         Min_Diff = 1 - Min_Perc,
         Max_Diff = 1 - Max_Perc) %>%
  mutate_if(is.numeric, round, 2) %>%
  as.data.frame

# Kruskal-Wallis test for significant difference
dci.kw <- kruskal.test(DCI_p ~ Method, data = data.dci)

# Dunn test for pairwise significant difference
dci.dunn <- dunnTest(DCI_p ~ Method, data = data.dci,
                     method = 'bonferroni')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# FIGURES: Figure ## - Model accuracy and variable importance
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Set inputs
data.model <- import_list('Results_GBM_Models.xlsx')

# Separate datasets
data.accuracy <- data.model$Model_Stats %>% rename(Method = Model) %>%
                                            mutate(Method = ifelse(Method == 'Passability', 'Rating category',
                                                                   'Presence/absence'),
                                                   AUC.rnd = round(AUC, 2))
data.variables <- data.model$Variable_Importance%>% rename(Method = Model) %>%
  mutate(Method = ifelse(Method == 'Passability', 'Rating category',
                         'Presence/absence'))

# Rename predictor variables
data.variables <- data.variables %>% mutate(Variable_Name = ifelse(var == 'Elev_m', 'Elevation',
                                                            ifelse(var == 'Slope_site_perc', 'Site slope',
                                                            ifelse(var == 'Slope_reach_perc', 'Reach slope',
                                                            ifelse(var == 'Slope_segment_perc', 'Segment slope',
                                                            ifelse(var == 'Qmad_cfs', 'Discharge',
                                                            ifelse(var == 'tnmfrc', 'Road class', 
                                                                   NA))))))) 

# Calculate stats for variable importance
stat.var <- data.variables %>% group_by(Method, var) %>%
  summarize(n = n(),
            Mean = mean(rel.inf, na.rm = TRUE),
            Min = min(rel.inf),
            Max = max(rel.inf),
            Median = median(rel.inf, na.rm = TRUE),
            SD = sd(rel.inf, na.rm = TRUE)) %>% 
  mutate(SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE),
         Range = Max - Min) %>%
  mutate_if(is.numeric, round, 2) %>%
  arrange(desc(Mean), .by_group = TRUE) %>%
  as.data.frame

# Calculate stats for model accuracy
stat.acc <- data.accuracy %>% group_by(Method) %>%
  summarize(n = n(),
            Mean = mean(AUC.rnd, na.rm = TRUE),
            Median = median(AUC, na.rm = TRUE),
            SD = sd(AUC, na.rm = TRUE),
            Min = min(AUC),
            Max = max(AUC),
            n_gt50 = sum(AUC.rnd > 0.5),
            n_gt70 = sum(AUC.rnd >= 0.7)) %>% 
  mutate(SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE)) %>%
  mutate_if(is.numeric, round, 2) %>%
  as.data.frame

# Calculate where mean and 95%ci lines intersect curves
dens.acc <- ggplot_build(ggplot(data.accuracy, aes(x = AUC, colour = Method)) + geom_density(aes(y = after_stat(scaled))))$data[[1]] %>%
  # mutate(Method = ifelse(group == 1, 'Presence/absence',
  #                        ifelse(group == 2, 'Landscape predictor', 'Random'))) %>%
  mutate(Method = ifelse(group == 1, 'Presence/absence',
                         ifelse(group == 2, 'Rating category', 'Random'))) %>%
  left_join(stat.acc) %>%
  select(x, y, Method, Mean, ci95.lo, ci95.up) %>%
  group_by(Method) %>%
  mutate(dens.mean = approx(x, y, xout = Mean)[[2]],
         dens.cilo = approx(x, y, xout = ci95.lo)[[2]],
         dens.cihi = approx(x, y, xout = ci95.up)[[2]]) %>%
  select(-x, -y) %>%
  slice(1) %>%
  as.data.frame

# Set fonts
font_import(pattern = 'Times New Roman')
loadfonts(device = 'win')

# Set colors
methods <- c('Uniform', 'Presence/absence', 'Rating category', 'Random')
methods.color <- viridis(n = length(methods), option = 'viridis', begin = 0.15)
names(methods.color) <- methods

# Plot AUC scores
# windows()
plot.accuracy <-
ggplot(data = data.accuracy, aes(x = AUC)) +
  geom_density(aes(fill = Method, y = after_stat(scaled)),
               alpha = 0.4, color = 'black') +
  geom_segment(data = dens.acc, aes(x = Mean, xend = Mean, y = 0, yend = dens.mean),
               linetype = 'dotted', linewidth = 0.5, color = 'black') +
  # geom_rug(aes(color = Method),
  #          sides = 'top') +
  geom_text(x = 0.34, y = 0.955, aes(label = '(a)', family = "Times New Roman"),
            size = 10/.pt) +
  
  scale_color_manual('Method', values = methods.color, 
                     breaks = c('Presence/absence', 'Rating category')) +
  scale_fill_manual('Method', values = methods.color, 
                    breaks = c('Presence/absence', 'Rating category')) +
  
  coord_cartesian(clip = 'off') +
  
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, family = "Times New Roman",
                                 color = 'black'),
        axis.title = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_text(vjust = -3, family = "Times New Roman"),
        axis.title.y = element_text(vjust = 0, family = "Times New Roman"),
        # plot.margin = margin(0.1, 0.2, 0.1, 0.1, 'cm'),
        legend.position = 'none') +
  
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = c(seq(0,1, by = 0.25))) +
  scale_x_continuous(limits = c(0.3,0.9), expand = c(0,0), breaks = c(seq(0,1, by = 0.1))) +
  labs(x = 'AUROC', y = 'Density')

# Plot Variable importance
data.variables$Method <- factor(data.variables$Method, levels = c('Rating category', 'Presence/absence'))
# windows()
plot.varimp <-
ggplot(data = data.variables, aes(x = rel.inf, y = reorder(Variable_Name, rel.inf, mean), fill = Method)) +
  stat_boxplot(geom = 'errorbar', width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(outlier.shape = NA, position = 'dodge', width = 0.5) +
  geom_text(x = 2.85, y = 6.3, aes(label = '(b)', family = "Times New Roman"),
            size = 10/.pt) +
  
  scale_color_manual('Method', values = methods.color, 
                     breaks = c('Presence/absence', 'Rating category')) +
  scale_fill_manual('Method', values = methods.color, 
                    breaks = c('Presence/absence', 'Rating category')) +
  
  # coord_cartesian(clip = 'off') +
  
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA, linewidth = 1),
        axis.text = element_text(size = 10, family = "Times New Roman",
                                 color = 'black'),
        axis.title = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_text(vjust = -1, family = "Times New Roman"),
        # plot.margin = margin(0.1, 0.1, 0.1, 0.1, 'cm'),
        legend.position = 'none') +
  
  scale_x_continuous(limits = c(-0.1,45), expand = c(0,0), breaks = c(seq(0, 50, by = 10))) +
  labs(x = 'Relative influence', y = "")

# Combine plots into single layout
# windows()
plot.out <-
plot_grid(plot.accuracy, plot.varimp,
          rel_widths = c(1,1.25))

# Create legend plot
plot.legend <- get_legend(plot.accuracy + guides(fill = guide_legend(nrow = 1,
                                                                     title = 'Alternative'))
                          + theme(legend.position = 'bottom',
                                  legend.text = element_text(size = 10, family = "Times New Roman"),
                                  legend.title = element_text(size = 10, family = 'Times New Roman'))) 

# Add legend
plot.out.out <-
plot_grid(plot.out, NULL, plot.legend,
          ncol = 1,
          rel_heights = c(1, 0.03, 0.1))


# Save output
ggsave('Figure_#_Model_performance.png', plot.out.out, 
       width = 16, height = 8, units = 'cm', dpi = 600, bg = 'white')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# FIGURES: Figure ## - Connectivity - Natural, dam, and road crossing
# ---------------------------------------------------------------------

# Set inputs
data.plot <- import_list('Results_Connectivity.xlsx')

# Combine distribution datasets
data.dci <- rbind(data.plot$Natural, data.plot$Dam, data.plot$Uniform, data.plot$PresAbs, data.plot$PassRate, data.plot$Random)

# Rename model scenarios
data.dci <- data.dci %>% mutate(Method = ifelse(Scenario == 'PresAbs', 'Presence/absence',
                                                ifelse(Scenario == 'PassRate', 'Rating category',
                                                       ifelse(Scenario == 'Random', 'Random sample', 
                                                              ifelse(Scenario == 'Uniform', 'Uniform', Scenario)))),
                                DCI_p = DCI_symm * 100)

# Calculate stats
stat.dci <- data.dci %>% group_by(Method) %>%
  summarize(n = n(),
            Mean = mean(DCI_p, na.rm = TRUE),
            Median = median(DCI_p, na.rm = TRUE),
            SD = sd(DCI_p, na.rm = TRUE),
            p_ltMean = sum(DCI_p > Mean) / n) %>% 
  mutate(SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE)) %>%
  mutate_if(is.numeric, round, 2) %>%
  as.data.frame

# Plotting data
methods <- c('Natural', 'Dam', 'Rating category', 'Presence/absence', 'Random sample', 'Uniform')
plot.dci <- stat.dci %>% select(Method, Mean) %>%
                         rename(DCIp = Mean) %>%
                         mutate(Barriers = ifelse(Method == 'Natural', 'Natural',
                                                  ifelse(Method == 'Dam', 'Dam', 'Road')),
                                DCIpL10 = log10(DCIp + 1),
                                middle = log10(DCIp + 1) / 2,
                                Barriers = factor(Barriers, levels = c('Natural', 'Dam', 'Road')),
                                Method = ordered(Method, levels = methods)) %>%
  arrange(Method)


# Set fonts
font_import()
loadfonts(device = 'win')

# Set colors
methods <- c('Natural', 'Dam', 'Uniform', 'Presence/absence', 'Rating category', 'Random sample')
methods.color <- c('steelblue3', 'gray50', viridis(n = 4, option = 'viridis', begin = 0.15))
names(methods.color) <- methods

# Plot
# windows()
plot.connect <-
ggplot(data = plot.dci, aes(x = Barriers, fill = Method)) +
  # geom_bar(stat = 'identity', position = 'identity')
  geom_tile(aes(height = DCIpL10, y = middle), width = 0.8) +
  scale_fill_manual('Alternative', values = methods.color, breaks = c('Uniform', 'Random sample',
                                                                      'Presence/absence',
                                                                      'Rating category')) +
  coord_cartesian(clip = 'off') +
  
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_text(family = "Times New Roman"),
        axis.title.y = element_text(family = "Times New Roman"),
        legend.position = c(0.68, 0.8),
        legend.text = element_text(size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 10, family = "Times New Roman"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,2),
                     expand = c(0,0),
                     breaks = c(0, 0.301, 1, 2),
                     labels = c(0, 1, 10, 100)) + 
  annotation_logticks(sides = 'l', outside = FALSE) +
  labs(y = bquote(DCI[p]), x = 'Barrier type')

# Save output
ggsave('Figure_#_Barrier_Connectivity.png', plot.connect, 
       width = 8, height = 10, units = 'cm', dpi = 600, bg = 'white')

# ---------------------------------------------------------------------


# --------------------------------------------------------------------- #
# FIGURES: Figure ## - Connectivity - Density plots
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Set inputs
data.plot <- import_list('Results_Connectivity.xlsx')

# Combine distribution datasets
data.dci <- rbind(data.plot$Uniform, data.plot$PresAbs, data.plot$PassRate, data.plot$Random)

# Rename model scenarios
data.dci <- data.dci %>% mutate(Method = ifelse(Scenario == 'PresAbs', 'Presence/absence',
                                                ifelse(Scenario == 'PassRate', 'Rating category',
                                                       ifelse(Scenario == 'Random', 'Random', 
                                                              ifelse(Scenario == 'Uniform', 'Uniform', NA)))),
                                DCI_p = DCI_symm * 100)

# Calculate stats
stat.dci <- data.dci %>% group_by(Method) %>%
  summarize(n = n(),
            Mean = mean(DCI_p, na.rm = TRUE),
            Median = median(DCI_p, na.rm = TRUE),
            SD = sd(DCI_p, na.rm = TRUE), 
            Min = min(DCI_p),
            Max = max(DCI_p)) %>% 
  mutate(SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE)) %>%
  mutate_if(is.numeric, round, 2) %>%
  as.data.frame

# Calculate where mean and 95%ci lines intersect curves
dens.dci <- ggplot_build(ggplot(data.dci, aes(x = DCI_p, colour = Method)) + geom_density(aes(y = after_stat(scaled))))$data[[1]] %>%
  mutate(Method = ifelse(group == 1, 'Presence/absence',
                         ifelse(group == 2, 'Random', 
                                ifelse(group == 3, 'Rating category', 'Uniform')))) %>%
  left_join(stat.dci) %>%
  select(x, y, Method, Mean, ci95.lo, ci95.up) %>%
  group_by(Method) %>%
  mutate(dens.mean = approx(x, y, xout = Mean)[[2]],
         dens.cilo = approx(x, y, xout = ci95.lo)[[2]],
         dens.cihi = approx(x, y, xout = ci95.up)[[2]]) %>%
  select(-x, -y) %>%
  slice(1) %>%
  as.data.frame

# Set labels
plot.labels <- data.frame(Method = c('Uniform', 'Random', 'Presence/absence', 'Rating category'),
                          label = c('(a) Uniform', 
                                    '(b) Random sample', 
                                    '(c) Presence/absence', 
                                    '(d) Rating category'))

# Set fonts
font_import()
loadfonts(device = 'win')

# Set colors
methods <- c('Uniform', 'Presence/absence', 'Rating category', 'Random')
methods.color <- viridis(n = length(methods), option = 'viridis', begin = 0.15)
names(methods.color) <- methods

# Ribbon 
ribbon.dat <- ggplot_build(ggplot(data = data.dci, aes(x = DCI_p, colour = Method)) +
                             geom_density(aes(y = after_stat(scaled))))$data[[1]] %>%
  mutate(Method = ifelse(group == 1, 'Presence/absence', 
                         ifelse(group == 2, 'Random',
                                ifelse(group == 3, 'Rating category', 'Uniform')))) %>%
  left_join(dens.dci) %>%
  group_by(Method) %>%
  filter(x >= ci95.lo & x <= ci95.up) %>%
  select(Method, x, y)

ribbon <- rbind(data.frame(Method = c('Presence/absence', 'Random', 'Rating category', 'Uniform'), 
                           x = c(0.77, 0.92, 1.17, 0.36), y = c(0,0,0,0)),
                as.data.frame(ribbon.dat),
                data.frame(Method = c('Presence/absence', 'Random', 'Rating category', 'Uniform'), 
                           x = c(1.33, 0.94, 2.04, 0.73), y = c(0,0,0,0)))

# Plot
# windows()
plot.dci <-
  ggplot() +
  geom_polygon(data = ribbon, aes(x = x, y = y, fill = Method), alpha = 0.6) +  
  
  geom_density(data = data.dci, aes(x = DCI_p, y = after_stat(scaled),
                                    fill = Method),
               alpha = 0.3, linewidth = 0.5, color = 'black') + 
  
  geom_segment(data = dens.dci, aes(x = Mean, xend = Mean, y = 0, yend = dens.mean),
               linetype = 'dotted', linewidth = 0.5, color = 'black') +
  
  geom_text(x = 1.05, y = 0.92, aes(label = label, family = windowsFont("Times New Roman")),
            data = plot.labels,
            hjust = 0) +
  
  # geom_rug(data = data.dci, aes(x = DCI_p, color = Method), 
  #          sides = 'top', length = unit(0.05, 'npc')) +
  
  facet_wrap(~fct_relevel(Method, 'Uniform', 'Random', 'Presence/absence', 'Rating category'), 
             ncol = 1) +
  
  scale_color_manual('Method', values = methods.color, breaks = c('Uniform', 'Random', 'Presence/absence', 'Rating category')) +
  scale_fill_manual('Method', values = methods.color, breaks = c('Uniform', 'Random', 'Presence/absence', 'Rating category')) +
  
  coord_cartesian(clip = 'off') +
  
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_text(vjust = -1, family = "Times New Roman"),
        axis.title.y = element_text(vjust = 3, family = "Times New Roman"),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, 'cm'), 
        # legend.position = c(0.75, 0.89),
        # legend.text = element_text(size = 10, family = "Times New Roman"),
        # legend.title = element_text(size = 10, family = "Times New Roman"),
        # legend.margin=margin(c(0,0,0,0)),
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        panel.spacing = unit(14, 'pt')) +
  
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = c(seq(0,1, by = 0.25))) +
  scale_x_continuous(limits = c(0,12), expand = c(0,0), breaks = c(seq(0,12, by = 1))) +
  labs(x = bquote(DCI[p]), y = 'Density')

# Save output
ggsave('Figure_#_Connectivity.png', plot.dci, width = 8.5, height = 16, units = 'cm', dpi = 600, bg = 'white')
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# FIGURES: Figure ## - Partial dependence plots
# ---------------------------------------------------------------------

# Load data
data.all <- rio::import_list('Results_GBM_Models.xlsx')

# Declare data, convert to natural units, and convert cfs to cms
data.pdp <- data.all[['PartialDependence']] %>%
  mutate(Value_nat = (10^Value) - 1) %>%
  mutate(Value_nat = ifelse(Variable == 'Qmad_cfs', Value_nat * 0.02832, Value_nat),
         Variable = ifelse(Variable == 'Qmad_cfs', 'Qmad_cms', Variable)) 

# Facet
windows()
ggplot(data = data.pdp, aes(x = Value_nat, y = yhat)) +
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0)) +
  facet_grid(Class ~ Variable, scales = 'free')

# Plots - Presence/absence
plot.elev <-
ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Elev_m'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Qmad_cms'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Slope_site_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Slope_reach_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Slope_segment_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1)
title <- 
ggdraw() +
  draw_label('(a) Presence', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.pa <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Impass
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Elev_m'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Qmad_cms'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Slope_site_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Slope_reach_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Slope_segment_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1)
title <- 
  ggdraw() +
  draw_label('(b) Impassable', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.impass <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Severe
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Elev_m'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Qmad_cms'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Slope_site_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Slope_reach_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Slope_segment_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1)
title <- 
  ggdraw() +
  draw_label('(c) Severe', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.severe <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Moderate
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Elev_m'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Qmad_cms'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Slope_site_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Slope_reach_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Slope_segment_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1)
title <- 
  ggdraw() +
  draw_label('(d) Moderate', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.moderate <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Passable
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Elev_m'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_text(size = 10, family = "Times New Roman")) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values', x = 'Elevation (m)')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Qmad_cms'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman")) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5)) +
  labs(x = bquote('Discharge ('*m^3~s^-1*')'))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Slope_site_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman")) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40)) +
    labs(x = 'Site slope (%)')

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Slope_reach_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman")) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10)) + 
  labs(x = 'Reach slope (%)')

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Slope_segment_perc'), aes(x = Value_nat, y = yhat))+
  geom_line(color = 'grey75', alpha = 0.5, aes(group = Model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman")) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10)) +
  labs(x = 'Segment slope (%)')


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1)
title <- 
  ggdraw() +
  draw_label('(e) Passable', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.passable <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Combine plots
plot.out <- plot_grid(plot.pa, plot.impass, plot.severe, plot.moderate, plot.passable,
                      ncol = 1)

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Write output
ggsave(filename = 'Figure_##_PartialDependence.png', plot.out,
       width = 18, height = 20, units = 'cm',
       dpi = 600, bg = 'white')

# ---------------------------------------------------------------------