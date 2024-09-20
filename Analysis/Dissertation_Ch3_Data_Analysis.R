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
if(!require("EnvStats")){
  install.packages("EnvStats", dependencies = TRUE); library(EnvStats)}
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
if(!require("PlaneGeometry")){
  install.packages("PlaneGeometry", dependencies = TRUE); library(PlaneGeometry)}
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
                                 field.ID){ # Field used to separate which records should be predicted
  
  # Initialize output list
  out.list <- list()
  
  # Generate seeds for reproducible random assignments
  set.seed(seed.start)
  seed.Random <- sample(x = 1000000, size = n, replace = FALSE)
  
  # Loop across all random assignment seeds
  for(i in 1:length(seed.Random)){
    # Initialize output data
    cols.out <- c(field.ID, field.Type, field.Pass)
    data.out <- data %>% select(one_of(cols.out)) %>%
      mutate(Pass_O = !!sym({{field.Pass}}))
    barriers.out <- data.frame(UID = character(),
                               Barrier_Expected = character(),
                               Pass = numeric(),
                               Pass_O = numeric(),
                               Pass_M = numeric(),
                               Pass_R = numeric())
    
    # Extract type field values to predict from 
    types.pred <- data.out %>% distinct(!!sym({{field.Type}})) %>% pull(!!sym({{field.Type}}))
    
    # Assign random passabilities by looping across each type
    for(j in 1:length(types.pred)){
      # Extract data by type
      data.type <- data.out %>% filter(!!sym({{field.Type}}) == types.pred[j])
      
      # Extract the passability ratings for a given type
      type <- freq.tbl %>% filter(!!sym({{field.Type}}) == types.pred[j]) %>% pull(!!sym({{field.Pass}}))
      
      # Extract passability rating probabilities for given type
      prob.type <- freq.tbl %>% filter(!!sym({{field.Type}}) == types.pred[j]) %>% pull(freq)
      
      # Randomly assign values to dataset
      set.seed(seed.Random[i])
      data.type <- data.type %>% mutate(Pass_M = sample(x = type,
                                                        size = nrow(data.type),
                                                        prob = prob.type,
                                                        replace = TRUE),
                                        Pass_R = ifelse(is.na(Pass), Pass_M, Pass)) %>%
        as.data.frame
      
      # Bind to output
      barriers.out <- rbind(barriers.out, data.type)
    }
    
    # Name output as seed
    name.seed <- paste0('Random_Seed_', as.character(seed.Random[i]))
    
    # Assign to output
    barriers.out <- barriers.out %>% mutate(Random_ID = name.seed)
    
    # Assign barrier randomization to output list
    out.list[[name.seed]] <- barriers.out
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
  
  # Model evaluation - TSS
  stat.cm <- as.data.frame(predictions.cm[['byClass']]) %>%
    mutate(Stat = rownames(.)) %>%
    pivot_wider(names_from = Stat, values_from = `predictions.cm[["byClass"]]`) %>% 
    mutate(TSS = Sensitivity + Specificity - 1) %>% as.data.frame()
  
  # Generate model evaluation dataframe
  tbl.eval <- data.frame(Model_seed = seed.split,
                         Accuracy = predictions.cm$overall[1],
                         Kappa = predictions.cm$overall[2],
                         AUC = predictions.auc$auc[1],
                         TSS = stat.cm$TSS,
                         row.names = NULL)
  
  # Generate names for list output
  name.model <- paste0('Model_Seed_', seed.split)
  name.eval  <- paste0('Eval_Seed_', seed.split)
  name.tss <- paste0('ClassStat_Seed_', seed.split)
  
  # Assign outputs to list
  out.list[[name.model]] <- model.out
  out.list[[name.eval]] <- tbl.eval
  out.list[[name.tss]] <- stat.cm
  
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
  
  # Model evaluation - TSS
  stat.cm <- as.data.frame(predictions.cm[['byClass']]) %>%
    mutate(Class = substring(rownames(.), 8),
           TSS = Sensitivity + Specificity - 1)
  row.names(stat.cm) <- NULL
  pred.tss <- mean(stat.cm[['TSS']], na.rm = T)
  
  # Generate model evaluation dataframe
  tbl.eval <- data.frame(Model_seed = seed.split,
                         Accuracy = predictions.cm$overall[1],
                         Kappa = predictions.cm$overall[2],
                         AUC = predictions.auc$auc[1],
                         TSS = pred.tss,
                         row.names = NULL)
  
  # Generate names for list output
  name.model <- paste0('Model_Seed_', seed.split)
  name.eval  <- paste0('Eval_Seed_', seed.split)
  name.tss <- paste0('ClassStat_Seed_', seed.split)
  
  # Assign outputs to list
  out.list[[name.model]] <- model.out
  out.list[[name.eval]] <- tbl.eval
  out.list[[name.tss]] <- stat.cm
  
  # Return output model
  return(out.list)
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
# Data: Initialize stream network and barrier data
# ---------------------------------------------------------------------

# Set inputs
file.data.stream <- 'Stream_BR.csv'

# Load in data
data.stream <- read.csv(file.data.stream,
                        header = TRUE)

# Declare barriers
data.barriers <- data.stream %>% filter(UID != "")

# Identify sampled barriers
data.sample <- data.stream %>% filter(Sampled == TRUE & Pass != '')

# Summarize barrier survey by sournce
source.sum <- data.barriers %>% group_by(SOURCE, Barrier_Expected) %>%
  summarise(n = n()) %>% as.data.frame

# Write output
data.all <- list()
data.all[['Streams']] <- data.stream
data.all[['Barriers']] <- data.barriers
data.all[['Barrier_Sample']] <- data.sample
data.all[['Source_Sum']] <- source.sum

# Export data
export(data.all, file = 'Data_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: Passability - Uniform Method
# ---------------------------------------------------------------------

# Summary: This section generates a barrier dataset with uniform passability ratings 

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- import_list(data.all.file)

# Initialize output list
barriers.uniform <- list()

# Define range of uniform passability values
pass.range <- c(0, 0.33, 0.67, 1)

# Uniformly assign passability ratings across range
for(i in 1:length(pass.range)){
  # Name scenario
  seed.id <- paste0('UniPass_', as.character(pass.range[i]))

  # Select barriers
  barriers.out <- data.barriers %>% select(UID, Barrier_Expected, Pass)
  
  # Add uniform passability rating and rating ID
  barriers.out <- barriers.out %>% mutate(Pass_0 = Pass,
                                          Pass_M = pass.range[i],
                                          Pass_R = ifelse((is.na(Pass) | Pass == ""),
                                                          pass.range[i],
                                                          Pass),
                                          Uniform_ID = seed.id)
  
  # Add to output list
  barriers.uniform[[seed.id]] <- barriers.out}

  # Combine to single dataset
  barriers.uniform.out <- do.call('rbind', barriers.uniform)
  
  # Calculate performance stats
  data.stat <- barriers.uniform.out %>%
    filter(Barrier_Expected != 'Waterfall' & Barrier_Expected != 'Dam' & !is.na(Pass))
  # Generate confusion matrix
  data.cm <- confusionMatrix(as.factor(data.stat$Pass_0), as.factor(data.stat$Pass_M))
  # Model evaluation - TSS
  stat.cm <- as.data.frame(data.cm[['byClass']]) %>%
    mutate(Class = substring(rownames(.), 8),
           TSS = Sensitivity + Specificity - 1,
           Method = 'Uniform',
           Model_seed = NA)
  # Aggregate for plotting
  tss.uniform <- stat.cm %>%
    group_by(Method, Model_seed) %>%
    dplyr::summarize(TSS = mean(TSS, na.rm = TRUE)) %>% as.data.frame
  
  # Write output
  data.all[['Barriers_Uniform']] <- barriers.uniform.out
  results[['Uniform_CStats']] <- stat.cm
  results[['Uniform_Stats']] <- tss.uniform
  
  # Export data
  export(data.all, file = 'Data_Connectivity.xlsx')
  export(results, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------
  
# --------------------------------------------------------------------- #
# METHODS: Passability - Random Sample
# ---------------------------------------------------------------------
  
  # Summary: This section generates the random type and passability assignment models.
  
  # Declare data
  data.all.file <- 'Data_Connectivity.xlsx'
  results.file <- 'Results_Connectivity.xlsx'
  
  # Load data
  data.all <- import_list(data.all.file)
  results <- import_list(results.file)
  data.sample <- data.all[['Barrier_Sample']]
  data.barriers <- data.all[['Barriers']]
  
  # Generate frequency table for expected barrier type
  freq.Type <- get_frequency_summary(data = data.sample, fields.Group = c('Barrier_Expected', 'Pass')) %>%
    add_row(Barrier_Expected = 'Waterfall', Pass = 0, n = 6, freq = 1.00) %>%
    add_row(Barrier_Expected = 'Dam', Pass = 0, n = 82, freq = 1.00) 
  
  
  # Generate random barrier type and passability rating assignments from observed distribution
  barriers.random <- generate_random_pass(data = data.barriers,
                                          freq.tbl = freq.Type,
                                          seed.start = 223445,
                                          n = 100,
                                          field.Type = 'Barrier_Expected',
                                          field.Pass = 'Pass',
                                          field.ID = 'UID')
  
  # Calculate performance stats for all random scenarios
  for(i in 1:length(barriers.random)){
    # Select sampled sites data
    loop.data <- barriers.random[[i]] %>%
      filter(Barrier_Expected != 'Waterfall' & Barrier_Expected != 'Dam' & !is.na(Pass))
    # Generate confusion matrix
    loop.cm <- confusionMatrix(as.factor(loop.data$Pass_O), as.factor(loop.data$Pass_M))
    # Model evaluation - TSS
    stat.cm <- as.data.frame(loop.cm[['byClass']]) %>%
      mutate(Class = substring(rownames(.), 8),
             TSS = Sensitivity + Specificity - 1,
             Method = 'Random',
             Model_seed = substring(names(barriers.random[i]), 13))
    row.names(stat.cm) <- NULL
    # Compile output
    if(i == 1){
      stat.out.all <- stat.cm
    } else {stat.out.all <- rbind(stat.out.all, stat.cm)}
  }
  
  # Calculate site TSS
  tss.site <- stat.out.all %>% 
    group_by(Method, Model_seed) %>%
    dplyr::summarize(TSS = mean(TSS, na.rm = TRUE)) %>% as.data.frame
  
  # Combine to single dataset
  barriers.random.out <- do.call('rbind', barriers.random)

  # Write output
  data.all[['Barriers_Random']] <- barriers.random.out
  results[['Random_CStats']] <- stat.out.all
  results[['Random_Stats']] <- tss.site

  # Export data
  export(data.all, file = 'Data_Connectivity.xlsx')
  export(results, file = 'Results_Connectivity.xlsx')
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: Passability - Presence/Absence machine learning method
# ---------------------------------------------------------------------

# Summary: This section generates the presence/absence barrier models.

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'
  
# Load data
data.all <- import_list(data.all.file)
results <- import_list(results.file)
data.sample <- data.all[['Barrier_Sample']]
data.barriers <- data.all[['Barriers']]
  
# Add barrier presence/absence field
data.sample.model <- data.sample %>% mutate(Pass_chr = paste0('P_', as.character(Pass)),
                                            Barrier_Present = ifelse(Pass_chr == 'P_1', 'Absent', 'Present'))

# Log10(n+1) transform continuous predictors
data.sample.log10 <- data.sample.model %>% mutate_at(.vars = c('Elev_m',
                                                                'Slope_site_perc',
                                                                'Slope_reach_perc',
                                                                'Slope_segment_perc',
                                                                'Qmad_cfs'),
                                                     .funs = function(x) log10(x + 1))

# Declare model inputs
# model.data.sample <-data.sample.log10
model.data.sample <-data.sample.model
model.response <- 'Barrier_Present'
model.predictors <- c('tnmfrc', 'Elev_m', 'Slope_site_perc',
                'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
model.perc.train <- 0.7
model.n.fold <- 5
model.type.sample <- 'down'
model.type.model <- 'gbm'
model.type.metric <- 'ROC'
model.seed.start <- 223445
model.n.models <- 100
model.type.observed  <- 'Barrier_Observed'
model.type.predicted <- 'Barrier_Expected'
model.pass.observed  <- 'Pass'
model.field.id       <- 'UID'

# Generate seeds
set.seed(model.seed.start)
seed.Model <- sample(x = 1000000, size = model.n.models, replace = FALSE)

# Generate blank list for models
models.out.all <- list()

# Generate model, evaluation stats, and predict barrier presence/absence.
for(i in 1:length(seed.Model)){
  # Generate binary model
  model.out <- generate_binary_model(data.sample = model.data.sample,
                                     response = model.response,
                                     predictors = model.predictors,
                                     seed.split = seed.Model[i],
                                     perc.train = model.perc.train,
                                     n.fold = model.n.fold,
                                     type.sample = model.type.sample,
                                     type.model = model.type.model,
                                     type.metric = model.type.metric)
  
  # Add model to output list
  models.out.all[[paste0('Binary_', as.character(seed.Model[i]))]] <- model.out[[1]]
  
  # Declare output data
  barriers.out <- data.barriers %>% select(all_of({{model.field.id}}),
                                           all_of({{model.type.observed}}),
                                           all_of({{model.type.predicted}}),
                                           all_of({{model.pass.observed}}),
                                           all_of({{model.predictors}}))
  
  # Initialize predicted passability fields
  barriers.out <- barriers.out %>% mutate(Pass_0 = Pass,
                                          Model_val = NA,
                                          Pass_M = NA,
                                          Pass_R = NA)
  
  # Assing passability ratings to dams and waterfalls
  barriers.out <- barriers.out %>% mutate(Pass_M = ifelse(Barrier_Expected == 'Dam' | Barrier_Expected == 'Waterfall',
                                                          0,NA),
                                          Pass_R = ifelse(Barrier_Expected == 'Dam' | Barrier_Expected == 'Waterfall',
                                                          0,NA))
  
  # Declare data for model prediction
  barriers.predict <- barriers.out %>% filter(is.na(Pass_M))
  
  # Predict barriers with model
  barriers.predict$Model_val <- predict(object = model.out[[1]],
                                        barriers.predict[,{{model.predictors}}],
                                        type = 'raw')
  
  # Assign passability based on barrier presence
  barriers.predict <- barriers.predict %>% mutate(Pass_M = ifelse(Model_val == 'Present', 0, 1),
                                                  Pass_R = ifelse(is.na(Pass), Pass_M, Pass_0))
  
  # Reduce data for joining
  barriers.join <- barriers.predict %>% select(all_of({{model.field.id}}), Model_val, Pass_M, Pass_R)
  
  # Join predictions to output barriers
  barriers.out <- barriers.out %>% left_join(barriers.join, by = model.field.id) %>%
    mutate(Model_val = if_else(is.na(Model_val.x), Model_val.y, Model_val.x),
           Pass_M = ifelse(is.na(Pass_M.x), Pass_M.y, Pass_M.x),
           Pass_R = ifelse(is.na(Pass_R.x), Pass_R.y, Pass_R.x),
           Binary_ID = paste0('Binary_Seed_', as.character(seed.Model[i]))) %>%
    select(-Model_val.x, -Model_val.y, -Pass_M.x, -Pass_M.y, -Pass_R.x, -Pass_R.y)
  
  # Add Seed ID to stats output
  model.out[[3]] <- model.out[[3]] %>% mutate(Model_seed = seed.Model[i])
  
  # Combine barrier prediction outputs
  if(i == 1){
    barriers.out.all <- barriers.out
  } else {barriers.out.all <- rbind(barriers.out.all, barriers.out)}
  
  if(i == 1){
    eval.out.all <- model.out[[2]]
  } else {eval.out.all <- rbind(eval.out.all, model.out[[2]])}
  
  if(i == 1){
    stat.out.all <- model.out[[3]]
  } else {stat.out.all <- rbind(stat.out.all, model.out[[3]])}
}
# Add method to eval
eval.out.all <- eval.out.all %>% mutate(Method = "Presence/Absence")
stat.out.all <- stat.out.all %>% mutate(Method = "Presence/Absence")

# Variable importance information for each model
varimp.list <- list()
for(i in 1:length(names(models.out.all))){
  data.out <- summary.gbm(models.out.all[[i]][[11]])
  varimp.list[[i]] <- data.out
}
eval.presabs.varimp <- as.data.frame(rbindlist(varimp.list, idcol = 'index')) %>%
  mutate(Method = 'Presence/Absence')

# Partial dependence plot (pdp) information from each model
pred.vars <- c('Elev_m', 'Qmad_cfs', 'Slope_site_perc', 'Slope_reach_perc', 'Slope_segment_perc')
names.presabs <- names(models.out.all)
# PDP - Initialize output
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
  for(j in 1:length(models.out.all)){
    # Retrieve data
    data.pdp <- pdp::partial(object = models.out.all[[j]],
                             pred.var = pred.vars[i],
                             which.class = 'Present',
                             plot = FALSE)
    # Format for output
    pdp.out <- data.pdp %>%
      rename_with(.cols = 1, ~'Value') %>%
      mutate(Model = names.presabs[j],
             Variable = pred.vars[i],
             Method = 'Presence/Absence',
             Class = 'Present')
    pdp.presabs <- rbind(pdp.presabs, pdp.out)
  }
}

# Write output
data.all[['Barriers_PresAbs']] <- barriers.out.all
results[['GBM_Stats_PresAbs']] <- eval.out.all
results[['GBM_CStats_PresAbs']] <- stat.out.all
results[['GBM_VarImp_PresAbs']] <- eval.presabs.varimp
results[['GBM_PDP_PresAbs']] <- pdp.presabs

# Export data
export(data.all, file = 'Data_Connectivity.xlsx')
export(results, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: Passability - Category machine learning method
# ---------------------------------------------------------------------

# Summary: This section generates the barrier passability rating models.

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
data.results.file <- 'Results_Connectivity.xlsx'
models.presabs <- list()

# Load data
data.all <- import_list(data.all.file)
results <- import_list(data.results.file)
data.sample <- data.all[['Barrier_Sample']]
data.barriers <- data.all[['Barriers']]

# Add barrier presence/absence field
data.sample.model <- data.sample %>% mutate(Pass_chr = paste0('P_', as.character(Pass)),
                                            Barrier_Present = ifelse(Pass_chr == 'P_1', 'Absent', 'Present'))

# Log10(n+1) transform continuous predictors
data.sample.log10 <- data.sample.model %>% mutate_at(.vars = c('Elev_m',
                                                               'Slope_site_perc',
                                                               'Slope_reach_perc',
                                                               'Slope_segment_perc',
                                                               'Qmad_cfs'),
                                                     .funs = function(x) log10(x + 1))

# Declare model inputs
# model.data.sample <-data.sample.log10
model.data.sample <-data.sample.model
model.response <- 'Pass_chr'
model.predictors <- c('tnmfrc', 'Elev_m', 'Slope_site_perc',
                      'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
model.perc.train <- 0.7
model.n.fold <- 5
model.type.sample <- 'up'
model.type.model <- 'gbm'
model.type.metric <- 'Accuracy'
model.seed.start <- 223445
model.n.models <- 100
model.type.observed  <- 'Barrier_Observed'
model.type.predicted <- 'Barrier_Expected'
model.pass.observed  <- 'Pass'
model.field.id       <- 'UID'

# Generate seeds
set.seed(model.seed.start)
seed.Model <- sample(x = 1000000, size = model.n.models, replace = FALSE)

# Generate blank list for models
models.out.all <- list()

# Generate model, evaluation stats, and predict barrier presence/absence.
for(i in 1:length(seed.Model)){
  # Generate binary model
  model.out <- generate_multiclass_model(data.sample = model.data.sample,
                                     response = model.response,
                                     predictors = model.predictors,
                                     seed.split = seed.Model[i],
                                     perc.train = model.perc.train,
                                     n.fold = model.n.fold,
                                     type.sample = model.type.sample,
                                     type.model = model.type.model,
                                     type.metric = model.type.metric)
  
  # Add model to output list
  models.out.all[[paste0('Binary_', as.character(seed.Model[i]))]] <- model.out[[1]]
  
  # Declare output data
  barriers.out <- data.barriers %>% select(all_of({{model.field.id}}),
                                           all_of({{model.type.observed}}),
                                           all_of({{model.type.predicted}}),
                                           all_of({{model.pass.observed}}),
                                           all_of({{model.predictors}}))
  
  # Initialize predicted passability fields
  barriers.out <- barriers.out %>% mutate(Pass_0 = Pass,
                                          Model_val = NA,
                                          Pass_M = NA,
                                          Pass_R = NA)
  
  # Assing passability ratings to dams and waterfalls
  barriers.out <- barriers.out %>% mutate(Pass_M = ifelse(Barrier_Expected == 'Dam' | Barrier_Expected == 'Waterfall',
                                                          0,NA),
                                          Pass_R = ifelse(Barrier_Expected == 'Dam' | Barrier_Expected == 'Waterfall',
                                                          0,NA))
  
  # Declare data for model prediction
  barriers.predict <- barriers.out %>% filter(is.na(Pass_M))
  
  # Predict barriers with model
  barriers.predict$Model_val <- predict(object = model.out[[1]],
                                        barriers.predict[,{{model.predictors}}],
                                        type = 'raw')
  
  # Assign passability based on barrier presence
  barriers.predict <- barriers.predict %>% mutate(Pass_M = ifelse(Model_val == 'P_0', 0, 
                                                                  ifelse(Model_val == 'P_0.33', 0.33,
                                                                         ifelse(Model_val == 'P_0.67', 0.67, 1))),
                                                  Pass_R = ifelse(is.na(Pass), Pass_M, Pass_0))
  
  # Reduce data for joining
  barriers.join <- barriers.predict %>% select(all_of({{model.field.id}}), Model_val, Pass_M, Pass_R)
  
  # Join predictions to output barriers
  barriers.out <- barriers.out %>% left_join(barriers.join, by = model.field.id) %>%
    mutate(Model_val = if_else(is.na(Model_val.x), Model_val.y, Model_val.x),
           Pass_M = ifelse(is.na(Pass_M.x), Pass_M.y, Pass_M.x),
           Pass_R = ifelse(is.na(Pass_R.x), Pass_R.y, Pass_R.x),
           Binary_ID = paste0('Category_Seed_', as.character(seed.Model[i]))) %>%
    select(-Model_val.x, -Model_val.y, -Pass_M.x, -Pass_M.y, -Pass_R.x, -Pass_R.y)
  
  # Add Seed ID to stats output
  model.out[[3]] <- model.out[[3]] %>% mutate(Model_seed = seed.Model[i])
  
  # Combine barrier prediction outputs
  if(i == 1){
    barriers.out.all <- barriers.out
  } else {barriers.out.all <- rbind(barriers.out.all, barriers.out)}
  
  if(i == 1){
    eval.out.all <- model.out[[2]]
  } else {eval.out.all <- rbind(eval.out.all, model.out[[2]])}
  
  if(i == 1){
    stat.out.all <- model.out[[3]]
  } else {stat.out.all <- rbind(stat.out.all, model.out[[3]])}
  
  # Print status
  print(paste0('Complete: ', as.character(seed.Model[i])))
}
# Add method to eval
eval.out.all <- eval.out.all %>% mutate(Method = "Category")
stat.out.all <- stat.out.all %>% mutate(Method = "Category")

# Variable importance information for each model
varimp.list <- list()
for(i in 1:length(names(models.out.all))){
  data.out <- summary.gbm(models.out.all[[i]][[11]])
  varimp.list[[i]] <- data.out
}
# Bind eval and stat lists into 
eval.category.varimp <- as.data.frame(rbindlist(varimp.list, idcol = 'index')) %>%
  mutate(Model = 'Category')

# Partial dependence plot (pdp) information from each model
pred.vars <- c('Elev_m', 'Qmad_cfs', 'Slope_site_perc', 'Slope_reach_perc', 'Slope_segment_perc')
names.category <- names(models.out.all)
class.pass <- c('P_0', 'P_0.33', 'P_0.67', 'P_1')
# PDP - Initialize output
data.out <- data.frame(Value = numeric(),
                       yhat = numeric(),
                       Model = character(),
                       Variable = character(),
                       Alternative = character(),
                       Class = character())
# Extract partial dependence plot info for each model and continuous predictor
## Presence/absence model
pdp.category <- data.out
for(i in 1:length(pred.vars)){
  for(j in 1:length(names.category)){
    for(k in 1:length(class.pass)){
      data.pdp <- pdp::partial(object = models.out.all[[j]],
                               pred.var = pred.vars[i],
                               which.class = class.pass[k],
                               plot = FALSE)
      # Format for output
      pdp.out <- data.pdp %>%
        rename_with(.cols = 1, ~'Value') %>%
        mutate(Model = names.category[j],
               Variable = pred.vars[i],
               Method = 'Category',
               Class = class.pass[k])
      pdp.category <- rbind(pdp.category, pdp.out)
    }
    # Print status
    print(paste0('Complete PDP: ', as.character(names.category[j]), ' ', pred.vars[i]))
  }
}

# Write output
data.all[['Barriers_Category']] <- barriers.out.all
results[['GBM_Stats_Category']] <- eval.out.all
results[['GBM_CStats_Category']] <- stat.out.all
results[['GBM_VarImp_Category']] <- eval.category.varimp
results[['GBM_PDP_Category']] <- pdp.category

# Export data
export(data.all, file = 'Data_Connectivity.xlsx')
export(results, file = 'Results_Connectivity.xlsx')
# --------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: Connectivity - Load network and generate igraph object
# NOTE: This step must be run prior to any connectivity analysis to initialize the graph object
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
#                          so DCI should equal 1 (i.e. no fragementation)
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
# METHODS: Connectivity - Natural barriers
# ---------------------------------------------------------------------

# Summary: Calculate DCI for natural barriers

# Read in results table as list
results.connectivity <- import_list('Results_Connectivity.xlsx')

# Set inputs
outScenario <- 'Connectivity_Natural'

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
# METHODS: Connectivity - Dams
# ---------------------------------------------------------------------

# Summary: Calculate DCI for dam barriers set to impassable (pass_u = pass_d = 0.0)

# Set inputs
outScenario <- 'Connectivity_Dam'

# Read in results table as list
results<- import_list('Results_Connectivity.xlsx')

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
results[[outScenario]] <- connectivity.out

# Save output
export(results, file = 'Results_Connectivity.xlsx')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# METHODS: Connectivity - Uniform method
# ---------------------------------------------------------------------

# Summary: Calculate DCI with uniform barrier passability rating assignments

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- import_list(data.all.file)
results.connectivity <- import_list(results.file)
barriers.uniform <- data.all[['Barriers_Uniform']]

# Set inputs
connectivity.barriers <- split(barriers.uniform, f = barriers.uniform$Uniform_ID) #Convert dataframe to list()
connectivity.scenario     <- 'Connectivity_Uniform'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

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
# METHODS: Connectivity - Random sample method
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- import_list(data.all.file)
barriers.random <- data.all[['Barriers_Random']]

# # Combine random barrier assemblages
# barriers.list <- list(barriers.random.1,barriers.random.2,barriers.random.3,barriers.random.4)
# barriers.random <- do.call('rbind', barriers.list)

# Set inputs
connectivity.barriers <- split(barriers.random, f = barriers.random$Random_ID) #Convert dataframe to list()
connectivity.scenario     <- 'Connectivity_Random'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

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
# METHODS: Connectivity - Presence/absence method
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- import_list(data.all.file)
results.connectivity <- import_list(results.file)
barriers.presabs <- data.all[['Barriers_PresAbs']]

# Set inputs
connectivity.barriers <- split(barriers.presabs, f = barriers.presabs$Binary_ID) #Convert dataframe to list()
connectivity.scenario     <- 'Connectivity_PresAbs'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

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
# METHODS: Connectivity - Rating category method
# ---------------------------------------------------------------------

# Summary: Calculate DCI for random barrier distributions

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- import_list(data.all.file)
results.connectivity <- import_list(results.file)
barriers.category <- data.all[['Barriers_Category']]

# Set inputs
connectivity.barriers <- split(barriers.category, f = barriers.category$Binary_ID) #Convert dataframe to list()
connectivity.scenario     <- 'Connectivity_Category'
connectivity.field.pass   <- 'Pass_R'
connectivity.field.weight <- 'Length_km'

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
# RESULTS: Barrier Surveys
# ---------------------------------------------------------------------

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'

# Load data
data.all <- import_list(data.all.file)
data.sample <- data.all[['Barrier_Sample']]

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
  mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame

# Summarize passability by limit
limit.sum <- data.sample %>% group_by(Limit_WDFG, Pass) %>%
  summarise(n = n()) %>%
  mutate(freq = round(n / sum(n), digits = 2)) %>% as.data.frame

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: Barrier method performance - TSS
# ---------------------------------------------------------------------

# Declare data
results.file <- 'Results_Connectivity.xlsx'

# Load data
results <- rio::import_list(results.file)
tss.uniform <- results[['Uniform_Stats']]
tss.random <- results[['Random_Stats']]
tss.presabs <- results[['GBM_Stats_PresAbs']]
tss.category <- results[['GBM_Stats_Category']]

# Combine TSS datasets
tss.uniform <- tss.uniform %>% mutate(Model_seed = as.character(Model_seed))
tss.presabs <- tss.presabs %>% select(Method, Model_seed, TSS) %>% mutate(Model_seed = as.character(Model_seed))
tss.category <- tss.category %>% select(Method, Model_seed, TSS) %>% mutate(Model_seed = as.character(Model_seed))
data.tss <- do.call('rbind', list(tss.uniform, tss.random, tss.presabs, tss.category))

# Set methods as ordered factors
data.tss <- data.tss %>%
  mutate(Method = ifelse(Method == 'Category', 'Rating category',
                         ifelse(Method == 'Presence/Absence', 'Presence/absence',
                                ifelse(Method == 'Random', 'Random sample', 'Uniform'))),
         Method = factor(Method, levels = c('Uniform', 'Random sample', 'Presence/absence', 'Rating category')))

# Calculate stats: TSS ~ Method
stat.tss <- data.tss %>%
  group_by(Method) %>%
  dplyr::summarize(mean = mean(TSS),
                   sd = sd(TSS),
                   min = min(TSS),
                   max = max(TSS),
                   q25 = quantile(TSS, p = 0.25),
                   q75 = quantile(TSS, p = 0.75),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         se.hi = mean + se,
         se.lo = mean -se,
         ci95 = qt(0.975, df = n-1) * se,
         ci95.lo = mean - ci95,
         ci95.hi = mean + ci95,
         pi95.lo = mean - (sd * 1.96),
         pi95.hi = mean + (sd * 1.96),
         range = max - min,
         range.pi = pi95.hi - pi95.lo) %>%
  replace(is.na(.), 0) %>%
  # mutate(across(where(is.numeric), round, 2)) %>%
  as.data.frame

# Plotting parameters
methods <- c('Uniform', 'Presence/absence', 'Rating category', 'Random sample')
methods.color <- viridis(n = length(methods), option = 'viridis', begin = 0.15)
names(methods.color) <- methods

# Plot TSS
plot.tss <-
  # windows()
  ggplot() +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_point(data = data.tss, aes(x = Method, y = TSS, group = Method, color = Method, fill = Method),
             shape = 21, alpha = 0.35, size = 0.75, position = position_jitter(w = 0.2, h = 0)) +
  geom_errorbar(data = stat.tss, aes(ymin = pi95.lo, ymax = pi95.hi, x = Method, group = Method), 
                width = 0.25) +
  geom_point(data = stat.tss, aes(y = mean, x = Method, fill = Method),
             shape = 21, color = 'black', size = 2.5) +
  scale_fill_manual(values = methods.color) +
  scale_color_manual(values = methods.color) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.x = element_text(size = 10, family = "Times New Roman", margin = unit(c(.25,0,0,0), 'cm')),
        axis.title.y = element_text(size = 10, family = "Times New Roman", margin = unit(c(0,.25,0,0), 'cm')),
        legend.position = 'none') +
  scale_y_continuous(limits = c(-0.25,0.65), expand = c(0,0), breaks = seq(-0.2,0.6,by=0.2)) +
  scale_x_discrete(labels = c('Uniform', 'Random\nsample', 'Presence/\nabsence', 'Rating\ncategory')) +
  labs(y = 'Predictive performance (TSS)')

# # Declare working directory
# pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
# setwd(pwd)
# 
# # Write output
# ggsave(filename = 'Figure_TSS.png', plot.tss,
#        width = 8.5, height = 10, units = 'cm',
#        dpi = 600, bg = 'white')
# 
# # Declare working directory
# pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: Barrier method performance - AUROC
# ---------------------------------------------------------------------

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- rio::import_list(data.all.file)
results <- rio::import_list(results.file)
barriers.uniform <- data.all[['Barriers_Uniform']] %>% mutate(Method = 'Uniform', 
                                                              Model_seed = readr::parse_number(Uniform_ID)) %>%
  rename(Pass_O = Pass_0) %>%
  select(-Uniform_ID)
barriers.random <- data.all[['Barriers_Random']] %>% mutate(Method = 'Random', 
                                                            Model_seed = readr::parse_number(Random_ID)) %>%
  select(-Random_ID)
stat.presabs <- results[['GBM_Stats_PresAbs']]
stat.category <- results[['GBM_Stats_Category']]

# AUROC - Uniform and Random Sample
barriers.all <- rbind(barriers.uniform, barriers.random)
barriers.sample <- barriers.all %>% filter(!is.na(Pass_O) & Barrier_Expected != 'Dam' & Barrier_Expected != 'Waterfall')
stat.barriers.sample <- data.frame(Model_seed = numeric(), Method = character(), AUC = numeric())
names.model <- unique(barriers.sample$Model_seed)
for(i in 1:length(names.model)){
  loop.data <- barriers.sample %>% filter(Model_seed == names.model[i])
  loop.auc <- multiclass.roc(loop.data$Pass_O, loop.data$Pass_M)
  loop.out <- data.frame(Model_seed = unique(loop.data$Model_seed),
                         Method = unique(loop.data$Method),
                         AUC = loop.auc$auc[1])
  stat.barriers.sample <- rbind(stat.barriers.sample, loop.out)
}

# Combine AUC datasets
data.auc <- do.call('rbind', list(stat.barriers.sample,
                                  stat.presabs %>% select(Model_seed, Method, AUC),
                                  stat.category %>% select(Model_seed, Method, AUC)))

# Set methods as ordered factors
data.auc <- data.auc %>%
  mutate(Method = ifelse(Method == 'Category', 'Rating category',
                         ifelse(Method == 'Presence/Absence', 'Presence/absence',
                                ifelse(Method == 'Random', 'Random sample', 'Uniform'))),
         Method = factor(Method, levels = c('Uniform', 'Random sample', 'Presence/absence', 'Rating category')))

# Calculate stats: TSS ~ Method
stat.auc <- data.auc %>%
  group_by(Method) %>%
  dplyr::summarize(mean = mean(AUC),
                   sd = sd(AUC),
                   min = min(AUC),
                   max = max(AUC),
                   q25 = quantile(AUC, p = 0.25),
                   q75 = quantile(AUC, p = 0.75),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         se.hi = mean + se,
         se.lo = mean -se,
         ci95 = qt(0.975, df = n-1) * se,
         ci95.lo = mean - ci95,
         ci95.hi = mean + ci95,
         pi95.lo = mean - (sd * 1.96),
         pi95.hi = mean + (sd * 1.96),
         range = max - min,
         range.pi = pi95.hi - pi95.lo) %>%
  # mutate(across(where(is.numeric), round, 2)) %>%
  as.data.frame

# Plotting parameters
methods <- c('Uniform', 'Presence/absence', 'Rating category', 'Random sample')
methods.color <- viridis(n = length(methods), option = 'viridis', begin = 0.15)
names(methods.color) <- methods

# Plot AUC
plot.auc <-
  #   windows()
  ggplot() +
  geom_hline(yintercept = 0.5, linetype = 'dotted') +
  geom_point(data = data.auc, aes(x = Method, y = AUC, group = Method, color = Method, fill = Method),
             shape = 21, alpha = 0.35, size = 0.75, position = position_jitter(w = 0.2, h = 0)) +
  geom_errorbar(data = stat.auc, aes(ymin = pi95.lo, ymax = pi95.hi, x = Method, group = Method),
                width = 0.25) +
  geom_point(data = stat.auc, aes(y = mean, x = Method, fill = Method),
             shape = 21, color = 'black', size = 2.5) +
  scale_fill_manual(values = methods.color) +
  scale_color_manual(values = methods.color) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.x = element_text(size = 10, family = "Times New Roman", margin = unit(c(.25,0,0,0), 'cm')),
        axis.title.y = element_text(size = 10, family = "Times New Roman", margin = unit(c(0,.25,0,0), 'cm')),
        legend.position = 'none') +
  scale_y_continuous(limits = c(0.38,0.82), expand = c(0,0), breaks = seq(0.4,0.8,by=0.1)) +
  scale_x_discrete(labels = c('Uniform', 'Random\nsample', 'Presence/\nabsence', 'Rating\ncategory')) +
  labs(y = 'Predictive performance (AUROC)')

# # Declare working directory
# pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
# setwd(pwd)
# 
# # Write output
# ggsave(filename = 'Figure_AUROC.png', plot.auc,
#        width = 8.5, height = 10, units = 'cm',
#        dpi = 600, bg = 'white')
# 
# # Declare working directory
# pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: Barrier method performance - TSS + AUROC
# ---------------------------------------------------------------------

# NOTE: Requires running above code for TSS and AUROC plots

# Combine TSS and AUROC plots
plot.performance <- 
# windows()
plot_grid(plot.tss, NULL, plot.auc, align = 'h',
          nrow = 1,
          rel_widths = c(1,0.1,1),
          labels = c('(a)', "", '(b)'),
          label_size = 10,
          label_fontfamily = 'Times New Roman',
          label_fontface = 'plain')

# Declare working directory
pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
setwd(pwd)

# Write output
ggsave(filename = 'Figure_4.pdf', plot.performance,
       device = cairo_pdf,
       width = 18, height = 10, units = 'cm',
       dpi = 800, bg = 'white')

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: GBM Performance - Partial Dependence Plots
# ---------------------------------------------------------------------

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
results.file <- 'Results_Connectivity.xlsx'

# Load data
data.all <- rio::import_list(data.all.file)
results <- rio::import_list(results.file)
data.pdp.presabs <- results[['GBM_PDP_PresAbs']]
data.pdp.category <- results[['GBM_PDP_Category']]
data.sample <- data.all[['Barrier_Sample']]
data.varimp <- rbind(results[["GBM_VarImp_PresAbs"]], results[["GBM_VarImp_Category"]])

# Average relative influence by Method
stat.varimp <- data.varimp %>%
  group_by(var, Method) %>%
  dplyr::summarize(mean = mean(rel.inf, na.rm = TRUE),
                   sd = sd(rel.inf),
                   min = min(rel.inf),
                   max = max(rel.inf),
                   n = n()) %>% 
  mutate(across(mean:max, round, 0)) %>%
  mutate(label = paste0('(',as.character(mean), '%)'))%>% as.data.frame()

# Rugs - as 10th percentiles 
data.sample.plot <- data.sample %>%
  select(UID,Elev_m, Slope_site_perc, Slope_reach_perc, Slope_segment_perc, Qmad_cfs, Pass) %>%
  mutate(Qmad_cms = Qmad_cfs * 0.02832) %>%
  pivot_longer(-c(UID, Pass), values_to = 'Value', names_to = 'Variable') %>%
  group_by(Variable) %>%
  dplyr::summarize(decile = scales::percent(c(seq(0.1,1,by = 0.1))),
                   Value = quantile(Value, seq(0.1,1,by=0.1))) %>%
  pivot_wider(names_from = 'Variable', values_from = 'Value') %>% as.data.frame
  
# Combine datasets
data.pdp <- rbind(data.pdp.presabs, data.pdp.category)

# Declare data, convert to natural units, and convert cfs to cms
data.pdp <- data.pdp %>%
  mutate(Value = ifelse(Variable == 'Qmad_cfs', Value * 0.02832, Value),
         Variable = ifelse(Variable == 'Qmad_cfs', 'Qmad_cms', Variable)) 

# Plotting parameters
rug.size <- 0.25
rug.length.cm <- 0.15
rug.side <- 't'
vars <- c('Elev_m', 'Qmad_cms', 'Slope_site_perc', 'Slope_reach_perc', 'Slope_segment_perc')
pa.model <- c('PA', 'Impassable', 'Severe', 'Moderate', 'Passable')
pa.breaks <- c('P', 'A')
pa.cols <- c('black', 'gray85')
line.col <- 'grey85'

# Plots - Presence/absence
plot.elev <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Elev_m'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  # geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf, color = !!sym({{pa.model[1]}})), 
  #          sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2250, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[1] & stat.varimp$Method == 'Presence/Absence')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Qmad_cms'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[2]}}), y = Inf), 
             sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2.9, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == 'Qmad_cfs' & stat.varimp$Method == 'Presence/Absence')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Slope_site_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[3]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 26, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[3] & stat.varimp$Method == 'Presence/Absence')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Slope_reach_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
    geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[4]}}), y = Inf), 
             sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[4] & stat.varimp$Method == 'Presence/Absence')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'Present' & Variable == 'Slope_segment_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[5]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[5] & stat.varimp$Method == 'Presence/Absence')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1, align = 'h', rel_widths = c(1.1,1,1,1,1))
title <- 
  ggdraw() +
  draw_label('(a) Presence', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.pa <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Impass
plot.elev <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Elev_m'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2250, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[1] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Qmad_cms'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[2]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2.9, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == 'Qmad_cfs' & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Slope_site_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[3]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 26, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[3] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Slope_reach_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[4]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[4] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0' & Variable == 'Slope_segment_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[5]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[5] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1, align = 'h', rel_widths = c(1.1,1,1,1,1))
title <- 
  ggdraw() +
  draw_label('(b) Impassable (0%)', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.impass <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Severe
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Elev_m'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2250, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[1] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Qmad_cms'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[2]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2.9, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == 'Qmad_cfs' & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Slope_site_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[3]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 26, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[3] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Slope_reach_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[4]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[4] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.33' & Variable == 'Slope_segment_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[5]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[5] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1, align = 'h', rel_widths = c(1.1,1,1,1,1))
title <- 
  ggdraw() +
  draw_label('(c) Severe (33%)', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.severe <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Moderate
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Elev_m'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2250, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[1] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Qmad_cms'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[2]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 2.9, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == 'Qmad_cfs' & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Slope_site_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[3]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 26, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[3] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Slope_reach_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[4]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[4] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment<-
  ggplot(data = data.pdp %>% filter(Class == 'P_0.67' & Variable == 'Slope_segment_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[5]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = 0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[5] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1, align = 'h', rel_widths = c(1.1,1,1,1,1))
title <- 
  ggdraw() +
  draw_label('(d) Moderate (67%)', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.moderate <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Plots - Passable
plot.elev <-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Elev_m'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = -0.55, x = 2250, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[1] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_text(size = 10, family = "Times New Roman"),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) +
  labs(y = 'fitted values', x = 'Elevation (m)')

plot.qmad <-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Qmad_cms'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[2]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = -0.55, x = 2.9, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == 'Qmad_cfs' & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman"),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5)) +
  labs(x = bquote('Discharge ('*m^3~s^-1*')'))

plot.sSite<-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Slope_site_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[3]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = -0.55, x = 26, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[3] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman"),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40)) +
  labs(x = 'Site slope (%)')

plot.sReach<-
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Slope_reach_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[4]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = -0.55, x = 7, hjust = 0, size = 10 / .pt,
           label = stat.varimp$label[which(stat.varimp$var == vars[4] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman"),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10)) + 
  labs(x = 'Reach slope (%)')

plot.sSegment<-
  # windows()
  ggplot(data = data.pdp %>% filter(Class == 'P_1' & Variable == 'Slope_segment_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[5]}}), y = Inf), 
           sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  annotate("text", y = -0.55, x = 7, hjust = 0, size = 10 / .pt,
            label = stat.varimp$label[which(stat.varimp$var == vars[5] & stat.varimp$Method == 'Category')],
           family = windowsFont("Times New Roman")) + 
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "Times New Roman"),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10)) +
  labs(x = 'Segment slope (%)')


plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1, align = 'h', rel_widths = c(1.1,1,1,1,1))
title <- 
  ggdraw() +
  draw_label('(e) Passable (100%)', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.passable <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

# Combine plots
plot.out <- plot_grid(plot.pa, plot.impass, plot.severe, plot.moderate, plot.passable,
                      ncol = 1)

# Declare working directory
pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
setwd(pwd)

# Write output
ggsave(filename = 'Figure_5.pdf', plot.out,
       device = cairo_pdf,
       width = 18, height = 20, units = 'cm',
       dpi = 600, bg = 'white')

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: Mean connectivity - Bar plot
# ---------------------------------------------------------------------

# Declare data
results.file <- 'Results_Connectivity.xlsx'

# Load data
results <- rio::import_list(results.file)

# Combine datasets
data.plot <- rbind(results$Connectivity_Natural, results$Connectivity_Dam, results$Connectivity_Uniform,
                   results$Connectivity_Random, results$Connectivity_PresAbs, results$Connectivity_Category)

# Rename model scenarios
data.dci <- data.plot %>% mutate(Method = ifelse(Scenario == 'Connectivity_PresAbs', 'Presence/absence',
                                                 ifelse(Scenario == 'Connectivity_Category', 'Rating category',
                                                        ifelse(Scenario == 'Connectivity_Random', 'Random sample', 
                                                               ifelse(Scenario == 'Connectivity_Uniform', 'Uniform', 
                                                                      ifelse(Scenario == 'Connectivity_Dam', 'Dams', 'Natural'))))),
                                 DCI_p = DCI_symm * 100) #%>%
#filter(Seed != 'UniPass_1' | is.na(Seed)) # Remove uniform passable scenario, which is same as dams

# Calculate stats
stat.dci <- data.dci %>% group_by(Method) %>%
  summarize(n = n(),
            Mean = mean(DCI_p, na.rm = TRUE),
            Median = median(DCI_p, na.rm = TRUE),
            SD = sd(DCI_p, na.rm = TRUE),
            Min = min(DCI_p),
            Max = max(DCI_p),
            qt.10 = quantile(DCI_p, probs = 0.1),
            qt.25 = quantile(DCI_p, probs = 0.25),
            qt.75 = quantile(DCI_p, probs = 0.75),
            p_ltMean = sum(DCI_p > Mean) / n) %>% 
  mutate(SE = SD / sqrt(n),
         ci95.lo = Mean - ((qt(0.975, df = n-1)) * SE),
         ci95.up = Mean + ((qt(0.975, df = n-1)) * SE),
         pi95 = SD * 1.96,
         pi95.lo = Mean - pi95,
         pi95.hi = Mean + pi95,
         Range = Max - Min) %>%
  mutate_if(is.numeric, round, 2) %>%
  as.data.frame

# Plotting data
methods <- c('Natural', 'Dams', 'Uniform', 'Rating category', 'Random sample', 'Presence/absence')
# methods <- c('Natural', 'Dam', 'Presence/absence', 'Rating category', 'Random sample', 'Uniform')
plot.dci <- stat.dci %>% select(Method, Mean, ci95.lo, ci95.up) %>%
  rename(DCIp = Mean) %>%
  mutate(Barriers = ifelse(Method == 'Natural', 'Natural',
                           ifelse(Method == 'Dams', 'Dams', 'Dams + Road')),
         DCIpL10 = log10(DCIp + 1),
         middle = log10(DCIp + 1) / 2,
         Barriers = factor(Barriers, levels = c('Natural', 'Dams', 'Dams + Road')),
         Method = ordered(Method, levels = methods),
         ci95.lo.l10 = log10(ci95.lo + 1),
         ci95.hi.l10 = log10(ci95.up + 1)) %>%
  mutate(ci95.lo.l10 = ifelse(is.na(ci95.lo.l10), 0.001, ci95.lo.l10)) %>%
  arrange(Method) 

# Set fonts
font_import()
loadfonts(device = 'win')

# Set colors
methods <- c('Natural', 'Dam', 'Uniform', 'Presence/absence', 'Rating category', 'Random sample')
methods.color <- c('steelblue3', 'gray50', viridis(n = 4, option = 'viridis', begin = 0.15))
names(methods.color) <- methods

seq.ticks <- data.frame(seq.nat = c(c(seq(0,1,by = 0.1)), c(seq(2,10,by = 1)), c(seq(20,100,by=10)))) %>%
  mutate(seq.l10 = log10(seq.nat + 1),
         label = ifelse(seq.nat == 0, 0,
                        ifelse(seq.nat == 1,1,
                               ifelse(seq.nat == 10, 10,
                                      ifelse(seq.nat == 100, 100, "")))))

# Plot
plot.out <-
  # windows()
  # ggplot(data = plot.dci, aes(x = Barriers, fill = Method)) +
  # # geom_bar(stat = 'identity', position = 'identity')
  # geom_tile(aes(height = DCIpL10, y = middle), width = 0.8) +
  ggplot(data = plot.dci, aes(y = DCIpL10, x = Barriers, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = ci95.lo.l10, ymax = ci95.hi.l10),
                position = position_dodge(0.9), 
                width = 0.4,
                linewidth = 0.25) +
  scale_fill_manual('Road-crossing barrier\nprediction method'
                    , values = methods.color, breaks = c('Uniform', 'Random sample',
                                                                            'Presence/absence',
                                                                            'Rating category')) +
  coord_cartesian(clip = 'off') +
  
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_text(family = "Times New Roman"),
        axis.title.y = element_text(family = "Times New Roman"),
        legend.position = c(0.63, 0.78),
        legend.text = element_text(size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 10, family = "Times New Roman", face = 'bold'),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,2.00432138), expand = c(0,0),
                     breaks = seq.ticks$seq.l10, label = seq.ticks$label) + 
  labs(x = 'Barrier type', y = bquote('Connectivity ('*DCI[p]*')'))

# Declare working directory
pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
setwd(pwd)

# Write output
ggsave(filename = 'Figure_6.pdf', plot.out,
       device = cairo_pdf,
       width = 8, height = 10, units = 'cm',
       dpi = 600, bg = 'white')

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: Method connectivity - Density plots
# ---------------------------------------------------------------------

# Declare data
results.file <- 'Results_Connectivity.xlsx'

# Load data
results <- rio::import_list(results.file)

# Combine datasets
data.plot <- rbind(results$Connectivity_Uniform, results$Connectivity_Random,
                   results$Connectivity_PresAbs, results$Connectivity_Category)

# Rename model scenarios
data.dci <- data.plot %>% mutate(Method = ifelse(Scenario == 'Connectivity_PresAbs', 'Presence/absence',
                                                 ifelse(Scenario == 'Connectivity_Category', 'Rating category',
                                                        ifelse(Scenario == 'Connectivity_Random', 'Random sample', 
                                                               ifelse(Scenario == 'Connectivity_Uniform', 'Uniform', 
                                                                      ifelse(Scenario == 'Connectivity_Dam', 'Dam', 'Natural'))))),
                                 DCI_p = DCI_symm * 100,
                                 DCI_l10 = log10(DCI_p + 1)) %>%
  filter(Seed != 'UniPass_1' | is.na(Seed)) # Remove uniform passable scenario, which is same as dams

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

# Log10 transform stats
stat.dci.l10 <- stat.dci %>%
  mutate(across(Mean:ci95.up, ~ log10(.x + 1)))

# Calculate where mean and 95%ci lines intersect curves
dens.dci <- ggplot_build(ggplot(data.dci, aes(x = DCI_l10, colour = Method)) + geom_density(aes(y = after_stat(scaled))))$data[[1]] %>%
  mutate(Method = ifelse(group == 1, 'Presence/absence',
                         ifelse(group == 2, 'Random sample',
                                ifelse(group == 3, 'Rating category', 'Uniform')))) %>%
  left_join(stat.dci.l10) %>%
  select(x, y, Method, Mean, ci95.lo, ci95.up) %>%
  group_by(Method) %>%
  mutate(dens.mean = approx(x, y, xout = Mean)[[2]],
         dens.cilo = approx(x, y, xout = ci95.lo)[[2]],
         dens.cihi = approx(x, y, xout = ci95.up)[[2]]) %>%
  select(-x, -y) %>%
  slice(1) %>%
  as.data.frame

# Generate ribbon for 95% CI
ribbon.dat <- ggplot_build(ggplot(data = data.dci, aes(x = DCI_l10, colour = Method)) +
                             geom_density(aes(y = after_stat(scaled))))$data[[1]] %>%
  mutate(Method = ifelse(group == 1, 'Presence/absence',
                         ifelse(group == 2, 'Random sample',
                                ifelse(group == 3, 'Rating category', 'Uniform')))) %>%
  left_join(dens.dci) %>%
  group_by(Method) %>%
  filter(x >= ci95.lo & x <= ci95.up) %>%
  select(Method, x, y)
ribbon <- rbind(data.frame(Method = c('Presence/absence', 'Random sample', 'Rating category', 'Uniform'),
                           x = c(0.2095150, 0.2810334, 0.4608978, 0.0211893), y = c(0,0,0,0)),
                as.data.frame(ribbon.dat),
                data.frame(Method = c('Presence/absence', 'Random sample', 'Rating category', 'Uniform'),
                           x = c(0.3304138, 0.2922561, 0.5943926, 0.1613680), y = c(0,0,0,0)))

# Set fonts
font_import()
loadfonts(device = 'win')

# Set colors
methods <- c('Uniform', 'Presence/absence', 'Rating category', 'Random')
methods.color <- viridis(n = length(methods), option = 'viridis', begin = 0.15)
names(methods.color) <- methods

seq.ticks <- data.frame(seq.nat = c(c(seq(0,1,by = 0.1)), c(seq(2,10,by = 1)), c(seq(20,100,by=10)))) %>%
  mutate(seq.l10 = log10(seq.nat + 1),
         label = ifelse(seq.nat == 0, 0,
                        ifelse(seq.nat == 1,1,
                               ifelse(seq.nat == 10, 10,
                                      ifelse(seq.nat == 100, 100, "")))))

# Plot - Uniform
plot.uniform <-
# windows()
  ggplot() +
  # geom_polygon(data = ribbon %>% filter(Method == 'Uniform'), aes(x = x, y = y),
  #              alpha = 0.6, fill = methods.color[1]) +
  geom_density(data = data.dci %>% filter(Method == 'Uniform'), aes(x = DCI_l10, y = after_stat(scaled)),
               alpha = 0.3, linewidth = 0.5, color = 'black', fill = methods.color[1]) + 
  geom_segment(data = dens.dci %>% filter(Method == 'Uniform'), aes(x = Mean, xend = Mean, y = 0, yend = dens.mean),
               linetype = 'dotted', linewidth = 0.5, color = 'black') +
  geom_text(y = 0.9, x = 0.4771213, hjust = 0, aes(label = '(a) Uniform', family = windowsFont("Times New Roman"))) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_blank(),
        # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = c(seq(0,1, by = 0.25))) +
  scale_x_continuous(limits = c(0,2.00432138), expand = c(0,0),
                     breaks = seq.ticks$seq.l10, label = seq.ticks$label) 

# Plot - Random sample
plot.random <-
# windows()
ggplot() +
  geom_polygon(data = ribbon %>% filter(Method == 'Random sample'), aes(x = x, y = y),
               alpha = 1, fill = methods.color[4]) +
  geom_density(data = data.dci %>% filter(Method == 'Random sample'), aes(x = DCI_l10, y = after_stat(scaled)),
               alpha = 0.3, linewidth = 0.5, color = 'black', fill = methods.color[4]) +
  geom_segment(data = dens.dci %>% filter(Method == 'Random sample'), aes(x = Mean, xend = Mean, y = 0, yend = dens.mean),
               linetype = 'dotted', linewidth = 0.5, color = 'black') +
  geom_text(y = 0.9, x = 0.4771213, hjust = 0, aes(label = '(b) Random sample', family = windowsFont("Times New Roman"))) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = c(seq(0,1, by = 0.25))) +
  scale_x_continuous(limits = c(0,2.00432138), expand = c(0,0),
                     breaks = seq.ticks$seq.l10, label = seq.ticks$label)

# Plot - Presence/absence
plot.presabs <-
# windows()
ggplot() +
  geom_polygon(data = ribbon %>% filter(Method == 'Presence/absence'), aes(x = x, y = y),
               alpha = 0.6, fill = methods.color[2]) +
  geom_density(data = data.dci %>% filter(Method == 'Presence/absence'), aes(x = DCI_l10, y = after_stat(scaled)),
               alpha = 0.3, linewidth = 0.5, color = 'black', fill = methods.color[2]) + 
  geom_segment(data = dens.dci %>% filter(Method == 'Presence/absence'), aes(x = Mean, xend = Mean, y = 0, yend = dens.mean),
               linetype = 'dotted', linewidth = 0.5, color = 'black') +
  geom_text(y = 0.9, x = 0.4771213, hjust = 0, aes(label = '(c) Presence/absence', family = windowsFont("Times New Roman"))) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = c(seq(0,1, by = 0.25))) +
  scale_x_continuous(limits = c(0,2.00432138), expand = c(0,0),
                     breaks = seq.ticks$seq.l10, label = seq.ticks$label) 

# Plot - Category
plot.category <-
# windows()
  ggplot() +
    geom_polygon(data = ribbon %>% filter(Method == 'Rating category'), aes(x = x, y = y),
                 alpha = 0.6, fill = methods.color[3]) +
    geom_density(data = data.dci %>% filter(Method == 'Rating category'), aes(x = DCI_l10, y = after_stat(scaled)),
                 alpha = 0.3, linewidth = 0.5, color = 'black', fill = methods.color[3]) + 
    geom_segment(data = dens.dci %>% filter(Method == 'Rating category'), aes(x = Mean, xend = Mean, y = 0, yend = dens.mean),
                 linetype = 'dotted', linewidth = 0.5, color = 'black') +
    geom_text(y = 0.9, x = 0.4771213, aes(label = '(d) Rating category', family = windowsFont("Times New Roman")), hjust = 0) +
    theme(panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(fill = NA),
          axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 0.5, 0, 0), 'cm'),
          legend.position = 'none') +
    coord_cartesian(clip = 'off') +
    scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = c(seq(0,1, by = 0.25))) +
  scale_x_continuous(limits = c(0,2.00432138), expand = c(0,0),
                     breaks = seq.ticks$seq.l10, label = seq.ticks$label) 

# Combine density plots to single figure
plot.dens <- plot_grid(plot.uniform, plot.random, plot.presabs, plot.category,
                       ncol = 1, align = 'v', rel_widths = c(1,1,1,1))

# Generate axis titles and blank space for final output
plot.xlab <- ggplot() + 
annotate('text', x = 0.5, y = 0.5, label = deparse(bquote('Connectivity ('*DCI[p]*')')), parse = TRUE,
         family = windowsFont("Times New Roman")) + theme_void()
plot.ylab <- ggplot() + 
  geom_text(y = 0.5, x = 0.5, angle = 90, aes(label = 'Density', family = windowsFont("Times New Roman"))) + theme_void()
plot.blank <- ggplot() + theme_void()

# Combine final plot
plot.out <- plot_grid(plot.ylab, plot.dens, plot.blank, plot.xlab,
                      nrow = 2, ncol = 2, 
                      rel_heights = c(1,0.05),
                      rel_widths = c(0.05,1))

# Declare working directory
pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
setwd(pwd)

# Write output
ggsave(filename = 'Figure_Connectivity_Density.pdf', plot.out,
       device = cairo_pdf,
       width = 8.5, height = 16, units = 'cm',
       dpi = 600, bg = 'white')

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# RESULTS: Method connectivity ~ TSS scatterplot
# ---------------------------------------------------------------------

# Declare data
results.file <- 'Results_Connectivity.xlsx'

# Load data
results <- rio::import_list(results.file)

# Combine connectivity datasets
data.conn <- rbind(results$Connectivity_Uniform, results$Connectivity_Random,
                   results$Connectivity_PresAbs, results$Connectivity_Category) %>%
  rename(Method = Scenario, 
         Model_seed = Seed) %>%
  mutate(Method = ifelse(Method == 'Connectivity_Category', 'Rating category',
                         ifelse(Method == 'Connectivity_PresAbs', 'Presence/absence',
                                ifelse(Method == 'Connectivity_Random', 'Random sample', 'Uniform'))),
         Model_seed = as.numeric(readr::parse_number(Model_seed)),
         Method = factor(Method, levels = c('Uniform', 'Random sample', 'Presence/absence', 'Rating category')),
         DCIp = DCI_symm * 100) %>%
  select(Method, Model_seed, DCIp)

# Combine TSS datasets
tss.uniform <- results[['Uniform_Stats']]
tss.random <- results[['Random_Stats']]
tss.presabs <- results[['GBM_Stats_PresAbs']]
tss.category <- results[['GBM_Stats_Category']]
tss.uniform <- tss.uniform %>% mutate(Model_seed = as.character(Model_seed))
tss.presabs <- tss.presabs %>% select(Method, Model_seed, TSS) %>% mutate(Model_seed = as.character(Model_seed))
tss.category <- tss.category %>% select(Method, Model_seed, TSS) %>% mutate(Model_seed = as.character(Model_seed))
data.tss <- do.call('rbind', list(tss.uniform, tss.random, tss.presabs, tss.category)) %>%
  filter(Method != 'Uniform') %>% 
  add_row(Method = c(rep('Uniform', 4)), Model_seed = c('0', '0.33', '0.67', '1'), TSS = c(rep(0,4))) %>%
  mutate(Method = ifelse(Method == 'Category', 'Rating category',
                         ifelse(Method == 'Presence/Absence', 'Presence/absence',
                                ifelse(Method == 'Random', 'Random sample', 'Uniform'))),
         Model_seed = as.numeric(readr::parse_number(Model_seed)),
         Method = factor(Method, levels = c('Uniform', 'Random sample', 'Presence/absence', 'Rating category')))

# Join connectivity and plotting datasets
data.plot <- data.conn %>%
  left_join(data.tss, by = c('Method', 'Model_seed')) %>%
  mutate(DCIp_Log10 = log10(DCIp + 1)) 

# Stats: TSS
stat.tss <- data.plot %>%
  filter(Model_seed != 1) %>%
  group_by(Method) %>%
  dplyr::summarize(mean = mean(TSS),
                   sd = sd(TSS),
                   min = min(TSS),
                   max = max(TSS),
                   q25 = quantile(TSS, p = 0.25),
                   q75 = quantile(TSS, p = 0.75),
                   q95.lo = quantile(TSS, p = 0.025),
                   q95.hi = quantile(TSS, p = 0.975),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci95 = qt(0.975, df = n-1) * se,
         ci95.lo = mean - ci95,
         ci95.hi = mean + ci95,
         pi95 = sd * 1.96,
         pi95.lo = mean - pi95,
         pi95.hi = mean + pi95,
         range = max - min,
         Metric = 'TSS') %>%
  replace(is.na(.), 0) %>%
  # mutate(across(where(is.numeric), round, 2)) %>%
  as.data.frame

# Stats: DCI
stat.dci <- data.plot %>%
  # filter(Model_seed != 1) %>%
  group_by(Method) %>%
  dplyr::summarize(mean = mean(DCIp),
                   sd = sd(DCIp),
                   min = min(DCIp),
                   max = max(DCIp),
                   q25 = quantile(DCIp, p = 0.25),
                   q75 = quantile(DCIp, p = 0.75),
                   q95.lo = quantile(DCIp, p = 0.025),
                   q95.hi = quantile(DCIp, p = 0.975),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci95 = qt(0.975, df = n-1) * se,
         ci95.lo = mean - ci95,
         ci95.hi = mean + ci95,
         pi95 = sd * 1.96,
         pi95.lo = mean - pi95,
         pi95.hi = mean + pi95,
         range = max - min,
         Metric = 'DCIp') %>%
  replace(is.na(.), 0) %>%
  # mutate(across(where(is.numeric), round, 2)) %>%
  as.data.frame

# Stats: TSS for plotting
stat.tss.plot <- stat.tss %>%
  select(Method, Metric, mean, ci95.lo, ci95.hi, pi95.hi)

# Stats: DCI for plotting
stat.dci.plot <- stat.dci %>%
  select(Method, Metric, mean, ci95.lo, ci95.hi, pi95.hi) %>%
  mutate(across(mean:pi95.hi, ~ log10(.x + 1)))

# Combine stats
stat.plot <- rbind(stat.tss.plot, stat.dci.plot) %>%
  pivot_longer(-c(Method, Metric), names_to = 'Stat', values_to = 'Value') %>%
  pivot_wider(names_from = c(Metric, Stat), values_from = Value) %>%
  mutate(TSS_pi95 = TSS_pi95.hi - TSS_mean,
         DCIp_pi95 = DCIp_pi95.hi - DCIp_mean) %>%
  as.data.frame %>%
  replace(is.na(.), 0.01)

# Ellipses
ell.random <- Ellipse$new(center = c(stat.plot$DCIp_mean[2], stat.plot$TSS_mean[2]),
                           rmajor = stat.plot$TSS_pi95[2], rminor = stat.plot$DCIp_pi95[2], alpha = 90)
ell.random.path <- ell.random $path()
ell.random.path <- as.data.frame(rbind(ell.random.path, ell.random.path[1,])) %>% 
  mutate(Method = 'Random sample')
ell.presabs <- Ellipse$new(center = c(stat.plot$DCIp_mean[3], stat.plot$TSS_mean[3]),
                           rmajor = stat.plot$DCIp_pi95[3], rminor = stat.plot$TSS_pi95[3], alpha = 0)
ell.presabs.path <- ell.presabs$path()
ell.presabs.path <- as.data.frame(rbind(ell.presabs.path, ell.presabs.path[1,])) %>% 
  mutate(Method = 'Presence/absence')
ell.category <- Ellipse$new(center = c(stat.plot$DCIp_mean[4], stat.plot$TSS_mean[4]),
                            rmajor = stat.plot$DCIp_pi95[4], rminor = stat.plot$TSS_pi95[4], alpha = 0)
ell.category.path <- ell.category$path()
ell.category.path <- as.data.frame(rbind(ell.category.path, ell.category.path[1,])) %>% 
  mutate(Method = 'Rating category')

data.ells <- do.call('rbind', list(ell.random.path,ell.presabs.path, ell.category.path))

# Set fonts
font_import()
loadfonts(device = 'win')

# Set colors
methods <- c('Uniform', 'Presence/absence', 'Rating category', 'Random sample')
methods.color <- viridis(n = length(methods), option = 'viridis', begin = 0.15)
names(methods.color) <- methods

seq.ticks <- data.frame(seq.nat = c(c(seq(0,1,by = 0.1)), c(seq(2,10,by = 1)), c(seq(20,100,by=10)))) %>%
  mutate(seq.l10 = log10(seq.nat + 1),
         label = ifelse(seq.nat == 0, 0,
                        ifelse(seq.nat == 1,1,
                               ifelse(seq.nat == 10, 10,
                                      ifelse(seq.nat == 100, 100, "")))))

# Plot
plot.out <-
# windows()
ggplot() +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_point(data = data.plot, 
             aes(y = TSS, x = DCIp_Log10, color = Method, fill = Method),
             alpha = 0.6, size = 0.6, shape = 21) +
  geom_path(data =  data.ells, aes(x = x, y = y, color = Method)) +
  geom_errorbar(data = stat.plot, 
                aes(xmin = DCIp_ci95.lo, xmax = DCIp_ci95.hi, y = TSS_mean, group = Method)) +
  geom_point(data = stat.plot, 
             aes(y = TSS_mean, x = DCIp_mean, fill = Method),
             alpha = 1, size = 2.75, shape = 21, color = 'black') +
  scale_fill_manual('Prediction method', values = methods.color) +
  scale_color_manual('Prediction method', values = methods.color) +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.x = element_text(size = 10, family = "Times New Roman", margin = unit(c(.25,0,0,0), 'cm')),
        axis.title.y = element_text(size = 10, family = "Times New Roman", margin = unit(c(0,.25,0,0), 'cm')),
        legend.position = c(0.72,0.83),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 10, family = "Times New Roman"),
        legend.title = element_text(size = 10, family = "Times New Roman", face = 'bold'),
        legend.key.height = unit(0.45, 'cm'),
        legend.key = element_blank(),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), 'cm')) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(-0.25,0.65), expand = c(0,0), breaks = seq(-0.2,0.6,by=0.2)) +
  scale_x_continuous(limits = c(0,2.00432138), expand = c(0,0),
                     breaks = seq.ticks$seq.l10, labels = seq.ticks$label) +
  labs(x = bquote('Connectivity ('*DCI[p]*')'), y = 'Predictive performance (TSS)') +
  guides(color = guide_legend(override.aes = list(linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0)))

# Declare working directory
pwd <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/Figures/')
setwd(pwd)

# Write output
ggsave(filename = 'Figure_7.pdf', plot.out,
       device = cairo_pdf,
       width = 8.5, height = 10, units = 'cm', 
       bg = 'white', dpi = 600)

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# APPENDIX: Rating category as continuous/regression ML model
# ---------------------------------------------------------------------

# NOTE: This analysis is used to evaluate multi-class partial dependence
#       plot trends by treating passability as a continuous class which
#       condenses partial dependence to singular plots. If plots match 
#       to categorical trends, relationships are likely a sound mechanistic
#       interpretation. If trends to not match categorical plots, then 
#       model type has influence and may require additional analysis.

# Declare data
data.all.file <- 'Data_Connectivity.xlsx'
data.results.file <- 'Results_Connectivity.xlsx'
models.presabs <- list()

# Load data
data.all <- import_list(data.all.file)
results <- import_list(data.results.file)
data.sample <- data.all[['Barrier_Sample']]
data.barriers <- data.all[['Barriers']]

# Add barrier presence/absence field
data.sample.model <- data.sample %>% mutate(Pass_chr = paste0('P_', as.character(Pass)),
                                            Barrier_Present = ifelse(Pass_chr == 'P_1', 'Absent', 'Present'))

# Declare model inputs
model.data.sample <-data.sample.model
model.response <- 'Pass'
model.predictors <- c('tnmfrc', 'Elev_m', 'Slope_site_perc',
                      'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
model.perc.train <- 0.7
model.n.fold <- 5
model.type.sample <- 'up'
model.type.model <- 'gbm'
model.type.metric <- 'Accuracy'
model.seed.start <- 223445
model.n.models <- 100
model.type.observed  <- 'Barrier_Observed'
model.type.predicted <- 'Barrier_Expected'
model.pass.observed  <- 'Pass'
model.field.id       <- 'UID'

# Generate seeds
set.seed(model.seed.start)
seed.Model <- sample(x = 1000000, size = model.n.models, replace = FALSE)

# Generate blank list for models
models.out.all <- list()

# Generate model, evaluation stats, and predict barrier presence/absence.
for(i in 1:length(seed.Model)){
  # Generate model
  # Initialize output list
  out.list <- list()
  # Select data for model
  data.model <- model.data.sample %>% select(all_of({{model.response}}), all_of({{model.predictors}})) %>%
    mutate_if(is.character, as.factor)
  # Split train/test datasets
  set.seed(seed.Model[i])
  train.Index <- createDataPartition(data.model[,{{model.response}}],
                                     p = model.perc.train,
                                     list = FALSE, times = 1)
  # Generate train/test split
  data.train <- data.model[train.Index,]
  data.test <- data.model[-train.Index,]
  
  # Tune model parameters
  # NOTE: Down/up-sampling is not applicable to regression models
  ctrl <- trainControl(method = 'cv', number = model.n.fold)
  
  # Declare formula
  model.formula <- reformulate('.', response = model.response)

  # Train model
  model.out <- train(model.formula,
                     data = data.train,
                     method = model.type.model,
                     trControl = ctrl,
                     metric = 'RMSE',
                     verbose = FALSE)
  
  models.out.all[[paste0('Model_Seed_', as.character(seed.Model[i]))]] <- model.out
}

# Partial dependence plot (pdp) information from each model
# pdp::partial(model.out, pred.var = 'Slope_site_perc', plot = TRUE, rug = TRUE)
pred.vars <- c('Elev_m', 'Qmad_cfs', 'Slope_site_perc', 'Slope_reach_perc', 'Slope_segment_perc')
names.category <- names(models.out.all)
# PDP - Initialize output
data.out <- data.frame(Value = numeric(),
                       yhat = numeric(),
                       Model = character(),
                       Variable = character())
# Extract partial dependence plot info for each model and continuous predictor
pdp.category <- data.out
for(i in 1:length(pred.vars)){
  for(j in 1:length(names.category)){
    data.pdp <- pdp::partial(object = models.out.all[[j]],
                             pred.var = pred.vars[i])
    # Format output
    pdp.out <- data.pdp %>%
      rename_with(.cols = 1, ~'Value') %>%
      mutate(Model = names.category[j],
             Variable = pred.vars[i])
    pdp.category <- rbind(pdp.category, pdp.out)
  }
  # Print status
  print(paste0('Complete PDP: ', as.character(names.category[j]), ' ', pred.vars[i]))
}

# Prep sample data for rugs
data.sample.plot <- data.sample %>%
  select(UID,Elev_m, Slope_site_perc, Slope_reach_perc, Slope_segment_perc, Qmad_cfs, Pass) %>%
  mutate(Qmad_cms = Qmad_cfs * 0.02832,
         PA = ifelse(Pass == 1, 'A', 'P'),
         Impassable = ifelse(Pass == 0, 'P', 'A'),
         Severe = ifelse(Pass == 0.33, 'P', 'A'),
         Moderate = ifelse(Pass == 0.67, 'P', 'A'),
         Passable = ifelse(Pass == 1, 'P', 'A')) %>%
  as.data.frame

# Declare data, convert to natural units, and convert cfs to cms
data.pdp <- pdp.category %>%
  mutate(Value = ifelse(Variable == 'Qmad_cfs', Value * 0.02832, Value),
         Variable = ifelse(Variable == 'Qmad_cfs', 'Qmad_cms', Variable))

# Plotting parameters
rug.size <- 0.15
rug.length.cm <- 0.15
rug.side <- 't'
vars <- c('Elev_m', 'Qmad_cms', 'Slope_site_perc', 'Slope_reach_perc', 'Slope_segment_perc')
pa.cols <- c('black', 'gray85')
line.col <- 'grey85'

# Plots - Presence/absence
plot.elev <-
  # windows()
  ggplot(data = data.pdp %>% filter(Variable == 'Elev_m'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  # geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[1]}}), y = Inf, color = !!sym({{pa.model[1]}})), 
  #          sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  # scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  # geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_text(size = 10, family = "Times New Roman"),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  # scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), breaks = c(1500, 2500)) 
  # labs(y = 'fitted values')

plot.qmad <-
  # windows()
  ggplot(data = data.pdp %>% filter(Variable == 'Qmad_cms'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  # geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[2]}}), y = Inf, color = !!sym({{pa.model[1]}})), 
  #          sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  # scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  # geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  # scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,5), breaks = c(0,2.5,5))

plot.sSite <-
  # windows()
  ggplot(data = data.pdp %>% filter(Variable == 'Slope_site_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  # geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[3]}}), y = Inf, color = !!sym({{pa.model[1]}})), 
  #          sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  # scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  # geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  # scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,45), breaks = c(0,20,40))

plot.sReach <-
  # windows()
  ggplot(data = data.pdp %>% filter(Variable == 'Slope_reach_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  # geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[4]}}), y = Inf, color = !!sym({{pa.model[1]}})),
  #          sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  # scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  # geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  # scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.sSegment <-
  # windows()
  ggplot(data = data.pdp %>% filter(Variable == 'Slope_segment_perc'), aes(x = Value, y = yhat))+
  geom_line(color = line.col, alpha = 0.5, aes(group = Model)) +
  # geom_rug(data = data.sample.plot, aes(x = !!sym({{vars[5]}}), y = Inf, color = !!sym({{pa.model[1]}})), 
  #          sides = rug.side, size = rug.size, length = unit(rug.length.cm, 'cm')) +
  # scale_color_manual(breaks = pa.breaks, values = pa.cols) +
  # geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
  geom_smooth(se = FALSE, color = 'black') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(size = 10, color = 'black', family = "Times New Roman"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  coord_cartesian(clip = 'off') +
  # scale_y_continuous(limits = c(-0.8, 0.8), expand = c(0,0), breaks = seq(-0.8, 0.8, by = 0.4)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = c(0,5,10))

plot.row <- plot_grid(plot.elev, plot.qmad, plot.sSite, plot.sReach, plot.sSegment,
                      nrow = 1, align = 'h', rel_widths = c(1.1,1,1,1,1))
title <- 
  ggdraw() +
  draw_label('(a) Regression', x = 0, hjust = 0, fontfamily = 'Times New Roman', size = 12) +
  theme(plot.margin = margin(0,0,0,7))

plot.regression <- plot_grid(title, plot.row, ncol = 1, rel_heights = c(0.15, 0.9))

windows()
plot.regression
# ---------------------------------------------------------------------