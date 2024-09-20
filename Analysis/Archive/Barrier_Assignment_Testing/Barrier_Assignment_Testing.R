# Title: Barrier_Assignment_Testing
# Author: Greg Goodrum
# Last update: 8/14/2023
# Contact: greg.goodrum@usu.edu
# Description: Testing workflow to calculate DCI from random assignment of barriers from observed data

# NOTES:
# Calculating connectivty indices
# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# References:
# https://cran.r-project.org/web/packages/riverconn/vignettes/Tutorial.html#generalized-riverscape-connectivity-index
# https://cran.r-project.org/web/packages/riverconn/riverconn.pdf
# https://damianobaldan.github.io/riverconn_tutorial/#reach-scale-indices
# https://igraph.org/r/doc/graph_from_data_frame.html
# https://igraph.org/r/doc/set_vertex_attr.html

# --------------------------------------------------------------------- #

# --------------------------------------------------------------------- #
# 01. Set up workspace
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

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 02. Initialize functions
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
# 03. Load data
# ---------------------------------------------------------------------

# Set inputs
raw.data            <- 'Stream_LR.csv'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv(raw.data,
                        header = TRUE)

# Assign passability of 0 for BHS when a water surface drop exists
data.stream <- data.stream %>% mutate(Pass_BHS = ifelse(!is.na(Water_surface_drop_m) & Water_surface_drop_m > 0,
                                                        '0.0', Pass_WDFG),
                                      Limit_BHS = ifelse(!is.na(Water_surface_drop_m) & Water_surface_drop_m > 0,
                                                         'WSD', Limit_WDFG))

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 04. Calculate barrier type and passability probabilities
# ---------------------------------------------------------------------

# Select sampled barriers
barrier.sample <- data.stream %>% filter(Sampled == 'TRUE' & Pass_WDFG != 'Unknown')

# Calculate probability of observed barrier types for each predicted type
prob.Type <- get_barrier_type_probability(data = barrier.sample,
                                          typePredicted = Barrier_Predicted,
                                          typeObserved = Barrier_Observed)

# Calculate probability of observed passability rating for each observed barrier type
prob.Pass <- get_barrier_pass_probability(data = barrier.sample,
                                          typeObserved = Barrier_Observed,
                                          passObserved = Pass_WDFG)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 05. Generate random barrier type and passability assignments based on observed frequency
# ---------------------------------------------------------------------

# Select barrier data and assign 'DAM' and 'FAL' an observed type
data.barriers <- data.stream %>% filter(grepl('UID', UID)) %>% 
                                 select(UID, Barrier_Predicted, Barrier_Observed, Pass_WDFG) %>% 
                                 mutate(Barrier_Observed = ifelse((Barrier_Predicted == 'DAM') | (Barrier_Predicted == 'FAL'),
                                        Barrier_Predicted,
                                        Barrier_Observed),
                                        Pass_WDFG = ifelse((Barrier_Predicted == 'DAM') | (Barrier_Predicted == 'FAL'),
                                                            '0.0',
                                                            Pass_WDFG))

# Generate random barrier type and passability data
barriers.random <- generate_random_barrier(data = data.barriers,
                                           seed = 1,
                                           nType = 3,
                                           nPass = 3,
                                           type_Pred = Barrier_Predicted,
                                           type_Obs = Barrier_Observed,
                                           type_Prob = prob.Type,
                                           pass_Obs = Pass_WDFG,
                                           pass_Prob = prob.Pass,
                                           type_Barr = Type_R,
                                           prob_Type = Barrier_Observed)

# Evaluate output type and pass frequencies
# data.check <- barriers.random$Type_Seed_234955$Pass_Seed_465902
# 
# check.Type_R <- get_barrier_type_probability(data = data.check,
#                                              typePredicted = Barrier_Predicted,
#                                              typeObserved = Type_R)
# 
# check.Pass_R <- get_barrier_pass_probability(data = data.check,
#                                              typeObserved = Type_R,
#                                              passObserved = Pass_R)

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 06. Create igraph network
# ---------------------------------------------------------------------

# Set inputs
graph.fields        <- c('From_Node', 'To_Node', 'Length_km', 'UID', 'Barrier_Predicted') # Fields to build graph with
attributes.edgeType <- 'Barrier_Predicted' # String corresponding to barrier type. 
attributes.edge     <- c('UID', 'Barrier_Predicted') # Strings corresponding to edge (barrier) attributes.
attributes.node     <- c('From_Node', 'Length_km') # String corresponding to node (stream segment) attributes. A field to weight connectivity indices is required, commonly assessed as length or HSI.
field.weight        <- 'Length_km' # String corresponding to node attribute used to weight connectivity indices.
file.graph          <- 'Graph_LR.png' # String indicating the name of the output file used to check igraph structure.

# Load in data
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

# Check graph network with dci calculation
# Initial output should return a dci = 1 where all barriers are passable.
# If the optional code is run, then the dci should be a numeric value between 0 and 1.
# OPTIONAL: Set barrier passability by barrier type
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   graph.stream <- set_edge_attr(graph.stream,
#                                 field.pass[i],
#                                 value = ifelse(E(graph.stream)$type == 'Confluence',
#                                                1.0,
#                                         ifelse(E(graph.stream)$type == 'Dam',
#                                                0.0,
#                                         ifelse(E(graph.stream)$type == 'Non-culvert Road Crossing',
#                                                0.75,
#                                         ifelse(E(graph.stream)$type == 'Unknown Road Crossing',
#                                                0.5,
#                                         ifelse(E(graph.stream)$type == 'Culvert',
#                                                0.25,
#                                                NA))))))
# }
index_calculation(graph = graph.stream,
                  weight = field.weight,
                  B_ij_flag = FALSE,
                  dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# ##. Join random barrier data to igraph
# ---------------------------------------------------------------------

# Calculate dci for all random barrier type and passability rating assignments
dci.random <- calculate_random_dci(inGraph = graph.stream,
                                   inBarriers = barriers.random,
                                   typeField = 'Type_R',
                                   passField = 'Pass_R',
                                   weightField = 'Length_km')

windows()
ggplot(dci.random, aes(x = DCI_symm)) +
  geom_density(fill = 'lightgray', aes(y = after_stat(scaled))) +
  theme_classic()

# ---------------------------------------------------------------------


View(as_data_frame(graph.stream, what = 'edges'))
View(as_data_frame(graph.stream, what = 'vertices'))

str(as_data_frame(graph.stream, what = 'edges'))
str(as_data_frame(graph.stream, what = 'vertices'))
