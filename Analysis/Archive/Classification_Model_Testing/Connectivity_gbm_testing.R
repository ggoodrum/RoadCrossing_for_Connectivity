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
# https://dpmartin42.github.io/posts/r/imbalanced-classes-part-1
# https://amunategui.github.io/binary-outcome-modeling/\


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
if(!require("lubridate")){
  install.packages("lubridate", dependencies = TRUE); library(lubridate)}
if(!require("igraph")){
  install.packages("igraph", dependencies = TRUE); library(igraph)}
if(!require("riverconn")){
  install.packages("riverconn", dependencies = TRUE); library(riverconn)}
if(!require("ggplot2")){
  install.packages("ggplot2", dependencies = TRUE); library(ggplot2)}
if(!require("ggnetwork")){
  install.packages("ggnetwork", dependencies = TRUE); library(ggnetwork)}
if(!require("viridis")){
  install.packages("viridis", dependencies = TRUE); library(viridis)}

if(!require("caret")){
  install.packages("caret", dependencies = TRUE); library(caret)}
if(!require("gbm")){
  install.packages("gbm", dependencies = TRUE); library(gbm)}
if(!require("pROC")){
  install.packages("pROC", dependencies = TRUE); library(pROC)}

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


# FUNCTION: get_barrier_type_probability() - Generate list of dataframes for each predicted barrier
#                                            type that compile the probability of a given observed
#                                            barrier type.

get_barrier_type_probability <- function(data,
                                         typePredicted,
                                         typeObserved){
  data.out <- data %>% group_by({{typePredicted}}, {{typeObserved}}) %>%
    summarise(n = n()) %>%
    mutate(freq = round(n /  sum(n), digits = 2))
  data.out <- as.data.frame(data.out)
  return(data.out)
}

# FUNCTION: get_barrier_pass_probability() - Generate list of dataframes for each observed barrier
#                                            type that compile the probability of a given observed
#                                            barrier passability rating.

get_barrier_pass_probability <- function(data,
                                         typeObserved,
                                         passObserved){
  data.out <- data %>% group_by({{typeObserved}}, {{passObserved}}) %>%
    summarise(n = n()) %>% 
    mutate(freq = round(n / sum(n), digits = 2))
  data.out <- as.data.frame(data.out)
  return(data.out)
}

# FUNCTION: generate_barrier_type() - Assign randomized barrier types to predicted instream barriers. Function
#                                     maintains observed barrier types (obsType), but assigns random type (Type_R)
#                                     to unsampled barriers based on observed type frequency (typeProb).

generate_barrier_type <- function(data, typePred, typeObs, typeProb, seed){
  # Initialize output dataframe
  data.out <- data %>% mutate(Type_R = NA)
  
  # Add observed barriers types to output
  data.out <- data.out %>% mutate(Type_R = ifelse(({{typeObs}} != "") & (!is.na({{typeObs}})),
                                                  {{typeObs}},
                                                  NA))
  
  # Generate random prediction dataset (sites lacking observation data)
  data.random <- data.out %>% filter(({{typeObs}} == "") | (is.na({{typeObs}})))
  
  # Identify all unobserved types
  types.Pred <- data.random %>% distinct({{typePred}}) %>% pull({{typePred}})
  
  # Generate seed list for random type samples
  set.seed(seed)
  seeds.Type <- sample(100000, length(types.Pred))
  
  # Loop across predicted barrier types and assign random type from observed frequency
  for(i in 1:length(types.Pred)){
    # Extract data for type
    data.type <- data.random %>% filter({{typePred}} == types.Pred[i])
    # Extract observed barrier types for random assignmnet
    type <- typeProb %>% filter({{typePred}} == types.Pred[i]) %>% pull({{typeObs}})
    # Extract observed barrier type probabilities for random assignment
    prob <- typeProb %>% filter({{typePred}} == types.Pred[i]) %>% pull(freq)
    # Randomly assign values
    set.seed(seeds.Type[i])
    data.type <- data.type %>% mutate(Type_R  = sample(x = type,
                                                       size = nrow(data.type),
                                                       prob = prob,
                                                       replace = TRUE)) %>%
      select(UID, Type_R)
    # Join to data.out
    data.out <- data.out %>% left_join(data.type, by = 'UID') %>%
      mutate(Type_R = coalesce(Type_R.y, Type_R.x)) %>%
      select(-Type_R.y, -Type_R.x)
  }
  
  # Return output data
  return(data.out)
}

# FUNCTION: generate_barrier_pass() - Assign randomized barrier passability rating for a given barrier type. Function
#                                     maintains observed passability rating(obsPass), but assigns a random
#                                     passability rating (Pass_R) to unsampled barriers based on a designated
#                                     barrier type (barrType) and observed passability rating frequency (passProb).

generate_barrier_pass <- function(data, typeBarr, passObs, passProb, probType, seed){
  # Initialize output dataframe
  data.out <- data %>% mutate(Pass_R = NA)
  
  # Add observed passability ratings to output
  data.out <- data.out %>% mutate(Pass_R = ifelse( ((({{passObs}} != "") & (!is.na({{passObs}}))) & ({{passObs}} != 'Unknown')),
                                                   {{passObs}},
                                                   NA))
  
  # Generate random prediction dataset (sites lacking observation data)
  data.random <- data.out %>% filter(is.na(Pass_R))
  
  # Identify all unobserved types
  types.Barr <- data.random %>% distinct({{typeBarr}}) %>% pull({{typeBarr}})
  
  # Generate seed list for random type samples
  set.seed(seed)
  seeds.Type <- sample(100000, length(types.Barr))
  
  # Loop across predicted barrier types and assign random passabilities from observed frequency
  for(i in 1:length(types.Barr)){
    # Extract data for type
    data.type <- data.random %>% filter({{typeBarr}} == types.Barr[i])
    # Extract observed barrier passability ratings for a given type for random assignment
    pass <- passProb %>% filter({{probType}} == types.Barr[i]) %>% pull({{passObs}})
    # Extract observed barrier passability rating probabilities for random assignment
    prob <- passProb %>% filter({{probType}} == types.Barr[i]) %>% pull(freq)
    # Randomly assign values
    set.seed(seeds.Type[i])
    data.type <- data.type %>% mutate(Pass_R  = sample(x = pass,
                                                       size = nrow(data.type),
                                                       prob = prob,
                                                       replace = TRUE)) %>%
      select(UID, Pass_R)
    # Join to data.out
    data.out <- data.out %>% left_join(data.type, by = 'UID') %>%
      mutate(Pass_R = coalesce(Pass_R.y, Pass_R.x)) %>%
      select(-Pass_R.y, -Pass_R.x)
  }
  
  # Return output data
  return(data.out)
}

# FUNCTION: generate_random_barrier() - Assigns repeatable random barrier types and passability ratings for unsampled
#                                       instream barriers from an observed distribution. Initial seed cascades and
#                                       determines all subsequent seed values, ensuring repeatability and reproducible
#                                       results.

generate_random_barrier <- function(data, 
                                    seed, # The starting seed from which all others will cascade
                                    nType, # Number of random type assignments to execute
                                    nPass, # Number of random passability rating assignments for each type assignment
                                    type_Pred, # Field name of the predicted barrier type
                                    type_Obs, # Field name of the observed barrier type
                                    type_Prob, # Data frame of observed barrier type frequencies
                                    pass_Obs, # Field name of observed passability rating
                                    pass_Prob, # Data frame of observed barrier passability rating frequencies
                                    type_Barr, # Field name of observed barrier type for passability assignment
                                    prob_Type){ # Field name of observed barrier type for passability assignment
  # Initialize output list
  out.list <- list()
  
  # Generate random seeds from which to draw type and pass seeds
  set.seed(seed)
  seed.start <- sample(x = 1000, size = 2, replace = FALSE)
  
  # Generate random type and pass seed vectors
  set.seed(seed.start[1])
  seed.Type <- sample(x = 1000000, size = nType, replace = FALSE)
  set.seed(seed.start[2])
  seed.Pass <- sample(x = 1000000, size = nPass, replace = FALSE)
  
  # Loop across all type and pass seeds to generate random barriers
  for(i in 1:length(seed.Type)){
    for(j in 1:length(seed.Pass)){
      # Generate random barrier types for unsampled barriers
      data.out.Type <- generate_barrier_type(data = data,
                                             typePred = {{type_Pred}},
                                             typeObs = {{type_Obs}},
                                             typeProb = type_Prob,
                                             seed = seed.Type[i])
      
      # Add random barrier passability ratings for unsampled barrier
      data.out.Type.Pass <- generate_barrier_pass(data = data.out.Type,
                                                  typeBarr = {{type_Barr}},
                                                  passObs = {{pass_Obs}},
                                                  passProb = pass_Prob,
                                                  probType = {{prob_Type}},
                                                  seed = seed.Pass[j])
      
      # Generate Type and Pass seed assignments for output
      name.Type <- paste0('Type_Seed_', as.character(seed.Type[i]))
      name.Seed <- paste0('Pass_Seed_', as.character(seed.Pass[j]))
      
      # Assign value to output list
      out.list[[name.Type]][[name.Seed]] <- data.out.Type.Pass
    }
  }
  
  # Return output data
  return(out.list)
} 

# FUNCTION: join_edge_attributes() - Function that joins edge attributes based on a common identifier field.
#                                    Function currently hard-coded for 'UID', 'Type_R', and 'Pass_R' fields.

join_edge_attributes <- function(inGraph, inData, typeField, passField){
  # Initialize output
  outGraph <- inGraph
  
  # Join type field
  outGraph <- set_edge_attr(outGraph,
                            name = as.name(typeField),
                            index = E(outGraph),
                            value = as.character(sapply(E(outGraph)$UID, function(x){
                              unlist(inData %>%
                                       filter(UID == x) %>%
                                       .[[as.name(typeField)]])
                            })))
  
  # Join pass field
  outGraph <- set_edge_attr(outGraph,
                            name = as.name(passField),
                            index = E(outGraph),
                            value = as.numeric(sapply(E(outGraph)$UID, function(x){
                              unlist(inData %>%
                                       filter(UID == x) %>%
                                       .[[as.name(passField)]])
                            })))
  
  # Convert pass field to numeric
  outGraph <- set_edge_attr(outGraph,
                            name = as.name(passField),
                            index = E(outGraph),
                            value = ifelse(E(outGraph)$type == 'Confluence',
                                           1.0, E(outGraph)$Pass_R))
  
  # Return output
  return(outGraph)
}

# FUNCTION: calculate_random_dci() - Function calculates symmetric and aysmmetric DCI from a provided igraph object,
#                                    list of dataframes of barrier passability, and fields indicating the type,
#                                    passability, and weight fields of the input tables.

calculate_random_dci <- function(inGraph, inBarriers, typeField, passField, weightField){
  # Initialize output table
  data.out <- data.frame(Seed_Type = character(0),
                         Seed_Pass = character(0),
                         DCI_symm = numeric(0),
                         DCI_asym = numeric(0))
  
  # Loop across all inBarriers
  for(i in 1:length(inBarriers)){
    for(j in 1:length(inBarriers[[i]])){
      # Initialize loop output
      loop.out <- data.frame(Seed_Type = character(0),
                             Seed_Pass = character(0),
                             DCI_symm = numeric(0),
                             DCI_asym = numeric(0)) 
      
      # Initialize loop graph
      graph.loop <- inGraph
      
      # Join single barrier edge
      graph.loop <- join_edge_attributes(inGraph = graph.loop,
                                         inData = inBarriers[[i]][[j]],
                                         typeField = 'Type_R',
                                         passField = 'Pass_R')
      
      # Convert to pass_u and pass_d fields
      field.pass <- c('pass_u', 'pass_d')
      for(k in 1:length(field.pass)){
        graph.loop <- set_edge_attr(graph.loop,
                                    field.pass[k],
                                    value = E(graph.loop)$Pass_R)
      }
      
      # Calculate symmetric DCI
      dci.symm <- index_calculation(graph = graph.loop,
                                    weight = weightField,
                                    B_ij_flag = FALSE,
                                    dir_fragmentation_type = 'symmetric')
      
      # Calculate asymmetric DCI
      dci.asym <- index_calculation(graph = graph.loop,
                                    weight = weightField,
                                    B_ij_flag = FALSE,
                                    dir_fragmentation_type = 'asymmetric')
      
      # Append data to loop output
      loop.out <- loop.out %>% add_row(Seed_Type =names(inBarriers)[[i]],
                                       Seed_Pass = names(inBarriers[[i]])[[j]],
                                       DCI_symm = dci.symm$index,
                                       DCI_asym = dci.asym$index)
      
      # Print current location
      print(paste0('COMPLETE: ', names(inBarriers)[[i]], ' (', i, ') ', 
                   ' ', names(inBarriers[[i]])[[j]], ' (', j, ')'))
      
      # Rbind loop to output
      data.out <- rbind(data.out, loop.out)
    }
  }
  
  # Return output
  return(data.out)
}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 01. Load and explore data.
# ---------------------------------------------------------------------

# Set inputs
raw.data  <- 'Stream_BR.csv'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv(raw.data,
                        header = TRUE)

# Generate barrier data
data.barriers <- data.stream %>% filter(UID != "")

# Extract sampled barriers
data.sample <- data.stream %>% filter(Sampled == TRUE & Pass_WDFG != 'Unknown')

# Summarized sampled barriers by presence/absence
data.sample <- data.sample %>% mutate(Barrier_Present = ifelse(Pass_WDFG == '1.0', 'Absent', 'Present'))

pres.sum <- data.sample %>% group_by(Barrier_Present) %>%
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

pres.sum

# Summarize sampled barrier types
pass.sum <- get_barrier_pass_probability(data = data.sample,
                                         # typeObserved = Barrier_Observed,
                                         passObserved = Pass_WDFG)
pass.sum

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 02a. Presence/absence model - caret pkg
# ---------------------------------------------------------------------

# Declare data
data.model <- data.sample %>% select(Barrier_Present, Barrier_Predicted,
                                     Elev_m, Slope_site_perc, Slope_reach_perc, Slope_segment_perc,
                                     Qmad_cfs) %>%
  mutate(Barrier_Predicted = as.factor(Barrier_Predicted),
         Barrier_Present = as.factor(Barrier_Present))

# Index to split 70%/30% train/test that maintains observed distribution of dependent variable
set.seed(451)
train.Index <- createDataPartition(data.model$Barrier_Present, p = 0.7,
                                   list = FALSE, times = 1)

# Generate train/test split
data.train <- data.model[train.Index,]
data.test <- data.model[-train.Index,]

# Tune model parameters
ctrl <- trainControl(method = 'cv', number = 5, returnResamp = 'none',
                     summaryFunction = twoClassSummary, classProbs = TRUE)

# Train model
model.gbm <- train(Barrier_Present ~ .,
                   data = data.train,
                   method = 'gbm',
                   trControl = ctrl,
                   metric = 'ROC')

# Test model - class value
predictions <- predict(object = model.gbm, data.test[,2:7], type = 'raw')
predictions.cm <- confusionMatrix(data.test$Barrier_Present, as.factor(predictions))
print(predictions.cm)

# Test model - probabilities
predictions.prob <- predict(object = model.gbm, data.test[,2:7], type = 'prob')
predictions.auc <- roc(ifelse(data.test[,"Barrier_Present"] == 'Present', 1, 0), predictions.prob[[2]])
print(predictions.auc$auc)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 02b. Presence/absence model + down/up-sampling - caret pkg
# ---------------------------------------------------------------------

# Declare data
data.model <- data.sample %>% select(Barrier_Present, Barrier_Predicted,
                                     Elev_m, Slope_site_perc, Slope_reach_perc, Slope_segment_perc,
                                     Qmad_cfs) %>%
  mutate(Barrier_Predicted = as.factor(Barrier_Predicted),
         Barrier_Present = as.factor(Barrier_Present))

# Index to split 70%/30% train/test that maintains observed distribution of dependent variable
set.seed(547)
# train.Index <- createDataPartition(data.model$Barrier_Present, p = 0.7,
#                                    list = FALSE, times = 1)
train.Index <- createDataPartition(data.model[,{{'Barrier_Present'}}], p = 0.7,
                                   list = FALSE, times = 1)

# Generate train/test split
data.train <- data.model[train.Index,]
data.test <- data.model[-train.Index,]

# Tune model parameters
ctrl <- trainControl(method = 'cv', number = 5, returnResamp = 'none',
                     summaryFunction = twoClassSummary, classProbs = TRUE)

# Generate training grid for hyperparameters
# NOTE: Grid tuning made no difference, and conforms to Han et al. (2021) which found hyperparameter
#       tuning to have no effect on model performance.
# caretGrid <- expand.grid(interaction.depth = c(1,3,5),
#                          n.trees = (0:10) * 10,
#                          shrinkage = c(0.1, 0.01, 0.001),
#                          n.minobsinnode = 10)

# Implement downsampling
ctrl$sampling <- 'down'

# Train model
model.gbm <- train(Barrier_Present ~ .,
                   data = data.train,
                   method = 'gbm',
                   trControl = ctrl,
                   # tuneGrid = caretGrid,
                   metric = 'ROC')

# Test model - class value
predictions <- predict(object = model.gbm, data.test[,2:7], type = 'raw')
predictions.cm <- confusionMatrix(data.test$Barrier_Present, as.factor(predictions))
print(predictions.cm)

# Test model - probabilities
predictions.prob <- predict(object = model.gbm, data.test[,2:7], type = 'prob')
predictions.auc <- roc(ifelse(data.test[,"Barrier_Present"] == 'Present', 1, 0), predictions.prob[[2]])
print(predictions.auc$auc)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 02c. Passability model + down/up-sampling - caret pkg
# ---------------------------------------------------------------------

# Declare data
data.model <- data.sample %>% select(Pass_WDFG, Barrier_Predicted,
                                     Elev_m, Slope_site_perc, Slope_reach_perc, Slope_segment_perc,
                                     Qmad_cfs) %>%
  mutate(Barrier_Predicted = as.factor(Barrier_Predicted),
         Pass_WDFG = as.factor(paste0('P_', Pass_WDFG)))

# Index to split 70%/30% train/test that maintains observed distribution of dependent variable
set.seed(547)
train.Index <- createDataPartition(data.model$Pass_WDFG, p = 0.7,
                                   list = FALSE, times = 1)

# Generate train/test split
data.train <- data.model[train.Index,]
data.test <- data.model[-train.Index,]

# Tune model parameters
ctrl <- trainControl(method = 'cv', number = 5, returnResamp = 'none',
                     summaryFunction = multiClassSummary, classProbs = TRUE)

# Generate training grid for hyperparameters
# NOTE: Grid tuning made no difference, and conforms to Han et al. (2021) which found hyperparameter
#       tuning to have no effect on model performance.
# caretGrid <- expand.grid(interaction.depth = c(1,3,5),
#                          n.trees = (0:10) * 10,
#                          shrinkage = c(0.1, 0.01, 0.001),
#                          n.minobsinnode = 10)

# Implement downsampling
ctrl$sampling <- 'up'

# Train model
model.gbm <- train(Pass_WDFG ~ .,
                   data = data.train,
                   method = 'gbm',
                   trControl = ctrl,
                   # tuneGrid = caretGrid,
                   metric = 'Accuracy')

# Test model - class value
predictions <- predict(object = model.gbm, data.test[,2:7], type = 'raw')
predictions.cm <- confusionMatrix(data.test$Pass_WDFG, as.factor(predictions))
print(predictions.cm)

# Test model - probabilities
predictions.prob <- predict(object = model.gbm, data.test[,2:7], type = 'prob')
predictions.auc <- multiclass.roc(data.test$Pass_WDFG, predictions.prob)
print(predictions.auc$auc)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 03a. FUNCTION: generate_binary_model
# ---------------------------------------------------------------------

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

# Test function
model.test <- generate_binary_model(data.sample = test.data,
                                    response = test.response,
                                    predictors = test.predictors,
                                    seed.split = test.seed.split,
                                    perc.train = test.perc.train,
                                    n.fold = test.n.fold,
                                    type.sample = test.type.sample,
                                    type.model = test.type.model,
                                    type.metric = test.type.metric)

# Set inputs
test.data <- data.sample
test.response <- 'Barrier_Present'
test.predictors <- c('Barrier_Predicted', 'Elev_m', 'Slope_site_perc',
                     'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
test.seed.split <- 917522
test.perc.train <- 0.7
test.n.fold <- 5
test.type.sample <- 'down'
test.type.model <- 'gbm'
test.type.metric <- 'ROC'

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 03b. FUNCTION: generate_multiClass_model
# ---------------------------------------------------------------------

generate_multiClass_model <- function(data.sample,
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

# Test function
model.test <- generate_multiClass_model(data.sample = test.data,
                                        response = test.response,
                                        predictors = test.predictors,
                                        seed.split = test.seed.split,
                                        perc.train = test.perc.train,
                                        n.fold = test.n.fold,
                                        type.sample = test.type.sample,
                                        type.model = test.type.model,
                                        type.metric = test.type.metric)

# Set inputs
test.data <- data.sample %>% mutate(Pass_WDFG = paste0('P_', Pass_WDFG))
test.response <- 'Pass_WDFG'
test.predictors <- c('Barrier_Predicted', 'Elev_m', 'Slope_site_perc',
                     'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
test.seed.split <- 917522
test.perc.train <- 0.7
test.n.fold <- 5
test.type.sample <- 'up'
test.type.model <- 'gbm'
test.type.metric <- 'ROC'

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 03a. FUNCTION: generate_modeled_barriers
# ---------------------------------------------------------------------

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
      
      # Assign 0.0 passability rating to dams and falls and maintain known passability ratings
      barriers.out <- barriers.out %>% mutate(Pass_M = ifelse((!!sym({{type.predicted}}) == 'DAM') | (!!sym({{type.predicted}}) == 'FAL'),
                                                               '0.0', 
                                                       ifelse((is.na(Pass_M)) & (!!sym({{pass.observed}}) != '') & (!is.na(!!sym({{pass.observed}}))) & (!!sym({{pass.observed}}) != 'Unknown'),
                                                              !!sym({{pass.observed}}), NA)))
      
      # Declare data for model prediction
      barriers.predict <- barriers.out %>% filter(is.na(Pass_M))
      
      # Predict barriers with model
      barriers.predict$Model_val <- predict(object = model.out[[1]],
                                            barriers.predict[,{{predictors}}],
                                            type = 'raw')
      
      # Assign passability based on barrier presence
      barriers.predict <- barriers.predict %>% mutate(Pass_M = ifelse(Model_val == 'Absent', '1.0', '0.5'))
      
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
      model.out <- generate_multiClass_model(data.sample = data.sample,
                                             response = response,
                                             predictors = predictors,
                                             seed.split = seed.Model[i],
                                             perc.train = perc.train,
                                             n.fold = n.fold,
                                             type.sample = type.sample,
                                             type.model = type.model,
                                             type.metric = type.metric)
      
      # Generate name assignments model and evaluation for output
      name.out <- paste0('multiclass_seed_', seed.Model[i])
      
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
      
      # Assign 0.0 passability rating to dams and falls and maintain known passability ratings
      barriers.out <- barriers.out %>% mutate(Pass_M = ifelse((!!sym({{type.predicted}}) == 'DAM') | (!!sym({{type.predicted}}) == 'FAL'),
                                                              '0.0', 
                                                              ifelse((is.na(Pass_M)) & (!!sym({{pass.observed}}) != '') & (!is.na(!!sym({{pass.observed}}))) & (!!sym({{pass.observed}}) != 'Unknown'),
                                                                     !!sym({{pass.observed}}), NA)))
      
      # Declare data for model prediction
      barriers.predict <- barriers.out %>% filter(is.na(Pass_M))
      
      # Predict barriers with model
      barriers.predict$Model_val <- predict(object = model.out[[1]],
                                            barriers.predict[,{{predictors}}],
                                            type = 'raw')
      
      # Assign passability based on barrier presence
      barriers.predict <- barriers.predict %>% mutate(Pass_M = as.character(gsub('P_', "", Model_val)))
      
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

# Test function
barriers.test <- generate_modeled_barriers(data.sample = test.data,
                                           response = test.response,
                                           predictors = test.predictors,
                                           perc.train = test.perc.train,
                                           n.fold = test.n.fold,
                                           type.sample = test.type.sample,
                                           type.model = test.type.model,
                                           type.metric = test.type.metric,
                                           
                                           data.predict = test.predict,
                                           seed.start = test.seed.start,
                                           n.models = test.n.models,
                                           str.model = test.str.model,
                                           
                                           type.observed = test.type.observed,
                                           type.predicted = test.type.predicted,
                                           pass.observed = test.pass.observed,
                                           field.id = test.field.id)

# # INPUTS: MULTICLASS MODEL TESTING
# test.data <- data.sample %>% mutate(Pass_WDFG = paste0('P_', Pass_WDFG))
# test.response <- 'Pass_WDFG'
# test.predictors <- c('Barrier_Predicted', 'Elev_m', 'Slope_site_perc',
#                      'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
# test.seed.split <- 687
# test.perc.train <- 0.7
# test.n.fold <- 5
# test.type.sample <- 'up'
# test.type.model <- 'gbm'
# test.type.metric <- 'ROC'
# 
# # Set inputs - generate_modeled_barriers
# test.predict <- data.barriers %>% select(UID, Barrier_Predicted, Barrier_Observed, Pass_WDFG,
#                                          Elev_m, Slope_site_perc, Slope_reach_perc, 
#                                          Slope_segment_perc, Qmad_cfs)
# test.seed.start <- 222
# test.n.models <- 2
# test.str.model <- 'multiclass'
# test.type.predicted <- 'Barrier_Predicted'
# test.type.observed <- 'Barrier_Observed'
# test.pass.observed <- 'Pass_WDFG'
# test.field.id <- 'UID'

# INPUTS: BINARY MODEL TESTING
test.data <- data.sample %>% mutate(Pass_WDFG = paste0('P_', Pass_WDFG))
test.response <- 'Barrier_Present'
test.predictors <- c('Barrier_Predicted', 'Elev_m', 'Slope_site_perc',
                     'Slope_reach_perc', 'Slope_segment_perc', 'Qmad_cfs')
test.seed.split <- 687
test.perc.train <- 0.7
test.n.fold <- 5
test.type.sample <- 'down'
test.type.model <- 'gbm'
test.type.metric <- 'ROC'

# Set inputs - generate_modeled_barriers
test.predict <- data.barriers %>% select(UID, Barrier_Predicted, Barrier_Observed, Pass_WDFG,
                                         Elev_m, Slope_site_perc, Slope_reach_perc,
                                         Slope_segment_perc, Qmad_cfs)
test.seed.start <- 222
test.n.models <- 2
test.str.model <- 'binary'
test.type.predicted <- 'Barrier_Predicted'
test.type.observed <- 'Barrier_Observed'
test.pass.observed <- 'Pass_WDFG'
test.field.id <- 'UID'

# ---------------------------------------------------------------------

