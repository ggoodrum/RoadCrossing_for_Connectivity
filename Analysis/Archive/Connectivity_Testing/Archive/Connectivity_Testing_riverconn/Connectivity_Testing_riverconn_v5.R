# Title: Connectivity_Testing
# Author: Greg Goodrum
# Last update: 7/14/2023
# Contact: greg.goodrum@usu.edu
# Description: Testing workflow for calculating stream network connectivity

# NOTES 
#   Terminal nodes: The input data processing workflow creates a non-data node at the outlet.
#                   riverconn uses reaches as nodes, so this node should not have a length and
#                   must be removed before creating the graph object. If set to 0, it's presence
#                   slightly increases the connectivity estimate.

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
if(!require("dplyr")){
  install.packages("dplyr"); library(dplyr)}
if(!require("igraph")){
  install.packages("igraph"); library(igraph)}
if(!require("lubridate")){
  install.packages("lubridate"); library(lubridate)}
if(!require("riverconn")){
  install.packages("riverconn"); library(riverconn)}
if(!require("ggplot2")){
  install.packages("ggplot2"); library(ggplot2)}
if(!require("ggnetwork")){
  install.packages("ggnetwork"); library(ggnetwork)}
if(!require("viridis")){
  install.packages("viridis"); library(viridis)}
if(!require("ggplot2")){
    install.packages("ggplot2"); library(ggplot2)}
if(!require("viridis")){
  install.packages("viridis"); library(viridis)}

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
  # data.edges <- data.edges %>% mutate(type = ifelse({{EdgeType_field}} == "",
  #                                                  'Confluence',
  #                                                  {{EdgeType_field}}))
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
         dpi = 1200)
  
  
  return(data.graph)
  
}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 03. Test Example
# ---------------------------------------------------------------------

# Set inputs
raw.data            <- 'Stream_Test.csv'
raw.fields          <- c('From_Node', 'To_Node', 'Shape_Leng', 'UID', 'BARRIER_GROUP', 'Passability')

attributes.edgeType <- 'BARRIER_GROUP'
attributes.edge     <- c('UID', 'BARRIER_GROUP', 'Passability')
attributes.node     <- c('From_Node', 'Shape_Leng')
file.graph          <- 'Graph_TEST.png'


# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv(raw.data,
                        header = TRUE) %>%
               select({{raw.fields}})

# Check network structure
check.data <- network_check(inData = data.stream,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node
outlet <- check.data$Terminal_Reaches$From_Node

# Generate igraph object
data.graph <- generate_attributed_igraph(inData          = data.stream,
                                         From_field      = From_Node,
                                         To_field        = To_Node,
                                         EdgeType_field  = {{attributes.edgeType}},
                                         Edge_attributes = {{attributes.edge}},
                                         Node_attributes = {{attributes.node}},
                                         Outlet_node     = outlet,
                                         graphFile       = file.graph)

# OPTIONAL: Set barrier passability by barrier type
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   data.graph <- set_edge_attr(data.graph,
#                                 field.pass[i],
#                                 value = ifelse(E(data.graph)$type == 'Confluence',
#                                                1.0,
#                                         ifelse(E(data.graph)$type == 'Dam',
#                                                0.0,
#                                         ifelse(E(data.graph)$type == 'Non-culvert Road Crossing',
#                                                0.75,
#                                         ifelse(E(data.graph)$type == 'Unknown Road Crossing',
#                                                0.5,
#                                         ifelse(E(data.graph)$type == 'Culvert',
#                                                0.25,
#                                                NA))))))
# }

# Calculate connectivty indices
# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity (symmetric DCI)
out.dci <- index_calculation(graph = data.graph,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 04. Logan River testing
# ---------------------------------------------------------------------

# Set inputs
raw.data            <- 'Stream_LR.csv'
raw.fields          <- c('From_Node', 'To_Node', 'Shape_Leng', 'UID', 'BARRIER_GROUP', 'Passability')

attributes.edgeType <- 'BARRIER_GROUP'
attributes.edge     <- c('UID', 'BARRIER_GROUP', 'Passability')
attributes.node     <- c('From_Node', 'Shape_Leng')
file.graph          <- 'Graph_LR.png'


# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv(raw.data,
                        header = TRUE) %>%
  select({{raw.fields}})

# Check network structure
check.data <- network_check(inData = data.stream,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node
outlet <- check.data$Terminal_Reaches$From_Node

# Generate igraph object
data.graph <- generate_attributed_igraph(inData          = data.stream,
                                         From_field      = From_Node,
                                         To_field        = To_Node,
                                         EdgeType_field  = {{attributes.edgeType}},
                                         Edge_attributes = {{attributes.edge}},
                                         Node_attributes = {{attributes.node}},
                                         Outlet_node     = outlet,
                                         graphFile       = file.graph)

# OPTIONAL: Set barrier passability by barrier type
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   data.graph <- set_edge_attr(data.graph,
#                                 field.pass[i],
#                                 value = ifelse(E(data.graph)$type == 'Confluence',
#                                                1.0,
#                                         ifelse(E(data.graph)$type == 'Dam',
#                                                0.0,
#                                         ifelse(E(data.graph)$type == 'Non-culvert Road Crossing',
#                                                0.75,
#                                         ifelse(E(data.graph)$type == 'Unknown Road Crossing',
#                                                0.5,
#                                         ifelse(E(data.graph)$type == 'Culvert',
#                                                0.25,
#                                                NA))))))
# }

# Calculate connectivty indices
# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity (symmetric DCI)
out.dci <- index_calculation(graph = data.graph,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 05. Upper Bear River testing
# ---------------------------------------------------------------------

# Set inputs
raw.data            <- 'Stream_UBR.csv'
raw.fields          <- c('From_Node', 'To_Node', 'Shape_Leng', 'UID', 'BARRIER_GROUP', 'Passability')

attributes.edgeType <- 'BARRIER_GROUP'
attributes.edge     <- c('UID', 'BARRIER_GROUP', 'Passability')
attributes.node     <- c('From_Node', 'Shape_Leng')
file.graph          <- 'Graph_UBR.png'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv(raw.data,
                        header = TRUE) %>%
  select({{raw.fields}})

# Check network structure
check.data <- network_check(inData = data.stream,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node
outlet <- check.data$Terminal_Reaches$From_Node

# Generate igraph object
data.graph <- generate_attributed_igraph(inData          = data.stream,
                                         From_field      = From_Node,
                                         To_field        = To_Node,
                                         EdgeType_field  = {{attributes.edgeType}},
                                         Edge_attributes = {{attributes.edge}},
                                         Node_attributes = {{attributes.node}},
                                         Outlet_node     = outlet,
                                         graphFile       = file.graph)

# OPTIONAL: Set barrier passability by barrier type
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   data.graph <- set_edge_attr(data.graph,
#                                 field.pass[i],
#                                 value = ifelse(E(data.graph)$type == 'Confluence',
#                                                1.0,
#                                         ifelse(E(data.graph)$type == 'Dam',
#                                                0.0,
#                                         ifelse(E(data.graph)$type == 'Non-culvert Road Crossing',
#                                                0.75,
#                                         ifelse(E(data.graph)$type == 'Unknown Road Crossing',
#                                                0.5,
#                                         ifelse(E(data.graph)$type == 'Culvert',
#                                                0.25,
#                                                NA))))))
# }

# Calculate connectivty indices
# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity (symmetric DCI)
out.dci <- index_calculation(graph = data.graph,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 06. Bear River testing
# ---------------------------------------------------------------------

# Set inputs
raw.data            <- 'Stream_BR.csv'
raw.fields          <- c('From_Node', 'To_Node', 'Shape_Leng', 'UID', 'BARRIER_GROUP', 'Passability')

attributes.edgeType <- 'BARRIER_GROUP'
attributes.edge     <- c('UID', 'BARRIER_GROUP', 'Passability')
attributes.node     <- c('From_Node', 'Shape_Leng')
file.graph          <- 'Graph_BR.png'

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv(raw.data,
                        header = TRUE) %>%
  select({{raw.fields}})

# Check network structure
check.data <- network_check(inData = data.stream,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node
outlet <- check.data$Terminal_Reaches$From_Node

# Generate igraph object
data.graph <- generate_attributed_igraph(inData          = data.stream,
                                         From_field      = From_Node,
                                         To_field        = To_Node,
                                         EdgeType_field  = {{attributes.edgeType}},
                                         Edge_attributes = {{attributes.edge}},
                                         Node_attributes = {{attributes.node}},
                                         Outlet_node     = outlet,
                                         graphFile       = file.graph)

# OPTIONAL: Set barrier passability by barrier type
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   data.graph <- set_edge_attr(data.graph,
#                               field.pass[i],
#                               value = ifelse(E(data.graph)$type == 'Confluence',
#                                              1.0,
#                                         ifelse(E(data.graph)$type == 'Dam',
#                                             0.0,
#                                         ifelse(E(data.graph)$type == 'Non-culvert Road Crossing',
#                                             0.75,
#                                         ifelse(E(data.graph)$type == 'Unknown Road Crossing',
#                                             0.5,
#                                         ifelse(E(data.graph)$type == 'Culvert',
#                                             0.25,
#                                         NA))))))
# }

# Calculate connectivty indices
# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity (symmetric DCI)
out.dci <- index_calculation(graph = data.graph,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

View(as_data_frame(data.graph, what = 'edges'))
View(as_data_frame(data.graph, what = 'vertices'))

str(as_data_frame(data.graph, what = 'edges'))
str(as_data_frame(data.graph, what = 'vertices'))
