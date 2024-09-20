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
  data.edges <- data.edges %>% mutate(type = ifelse({{EdgeType_field}} == "",
                                                   'Confluence',
                                                   {{EdgeType_field}}))
  
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
                                    data.nodes %>%
                                      filter(From_Node == x) %>%
                                      .[col]
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
  ifelse(length(V(data.graph)) <= 50, dim.Plot <- 10,
  ifelse(length(V(data.graph)) <= 100, dim.Plot <- 20,
  ifelse(length(V(data.graph)) <= 1000, dim.Plot <- 30,
  ifelse(length(V(data.graph)) <= 10000, dim.Plot <- 50,
         75))))
  
  # Plot to confirm
  gg0 <- ggnetwork(data.graph,
                   layout =  layout_as_tree(data.graph %>% as.undirected, root = index.outlet),
                   scale = FALSE)
  plot <- 
    ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_nodes(alpha = 0.3) +
    geom_edges(alpha = 0.5,
               arrow = arrow(length = unit(10, "pt"), type = "closed"),
               aes(color = type)) +
    scale_color_viridis(discrete = TRUE)+
    geom_nodetext(aes(label = name), fontface = "bold") +
    theme_blank()
  ggsave(graphFile, plot = plot, 
         width = dim.Plot, height = dim.Plot, units = 'in',
         dpi = 800)
  
  
  return(list(data.edges    = data.edges,
              data.nodes    = data.nodes,
              data.graph    = data.graph))
}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 03. Test Example
# ---------------------------------------------------------------------

# Clean workspace
# rm(list = ls())

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv('Stream_Test.csv',
                        header = TRUE) %>%
               select(From_Node, To_Node, Shape_Leng, UID, BARRIER_GROUP, Passability)

# Check network structure
check.data <- network_check(inData = data.stream,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node
outlet <- check.data$Terminal_Reaches$From_Node

# Generate igraph object
data.graph <- generate_attributed_igraph(inData = data.stream,
                                         From_field = From_Node,
                                         To_field = To_Node,
                                         EdgeType_field = BARRIER_GROUP,
                                         Edge_attributes = c(UID, BARRIER_GROUP, Passability),
                                         Node_attributes = c(From_Node, Shape_Leng),
                                         Outlet_node = outlet,
                                         graphFile = 'Graph_TEST.png')

# Extract graph object
graph.stream <- data.graph$data.graph

# Plot to confirm
gg0 <- ggnetwork(graph.stream,
                 layout =  layout_as_tree(graph.stream %>% as.undirected, root = 7),
                 scale = FALSE)

# windows()
plot <- 
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_nodes(alpha = 0.3) +
  geom_edges(alpha = 0.5,
             arrow = arrow(length = unit(10, "pt"), type = "closed"),
             aes(color = type)) +
  scale_color_viridis(discrete = TRUE)+
  geom_nodetext(aes(label = name), fontface = "bold") +
  theme_blank()
ggsave('Graph_TEST.png', plot = plot, 
       width = 10, height = 10, units = 'in',
       dpi = 800)


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

# Calculate connectivty indices

# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity
out.dci <- index_calculation(graph = graph.stream,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 04. UBR Test
# ---------------------------------------------------------------------

# Clean workspace
# rm(list = ls())

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv('Stream_UBR.csv',
                        header = TRUE) %>%
  select(From_Node, To_Node, Shape_Leng, UID, BARRIER_GROUP, Passability)

# Check network structure
check.data <- network_check(inData = data.stream,
                            From_field = From_Node,
                            To_field = To_Node)

# Set outlet node
outlet <- check.data$Terminal_Reaches$From_Node

# Generate igraph object
data.graph <- generate_attributed_igraph(inData = data.stream,
                                         From_field = From_Node,
                                         To_field = To_Node,
                                         EdgeType_field = BARRIER_GROUP,
                                         Edge_attributes = c(UID, BARRIER_GROUP, Passability),
                                         Node_attributes = c(Shape_Leng),
                                         Outlet_node = outlet,
                                         graphFile = 'Graph_UBR.png')

# Extract graph object
graph.stream <- data.graph$data.graph

# Plot to confirm
gg0 <- ggnetwork(graph.stream,
                 layout =  layout_as_tree(graph.stream %>% as.undirected, root = outlet),
                 scale = FALSE)

# windows()
plot <- 
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_nodes(alpha = 0.3) +
  geom_edges(alpha = 0.5,
             arrow = arrow(length = unit(10, "pt"), type = "closed"),
             aes(color = type)) +
  scale_color_viridis(discrete = TRUE)+
  geom_nodetext(aes(label = name), fontface = "bold") +
  theme_blank()
ggsave('Graph_TEST.png', plot = plot, 
       width = 30, height = 30, units = 'in',
       dpi = 800)

# Assign barrier passability  
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   graph.stream <- set_edge_attr(graph.stream,
#                                 field.pass[i],
#                                 value = ifelse(E(graph.stream)$type == 'Confluence',
#                                                1.0,
#                                                0.5))
# }

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

# Calculate connectivty indices

# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity
out.dci <- index_calculation(graph = graph.stream,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 05. LR Test
# ---------------------------------------------------------------------

# Clean workspace
# rm(list = ls())

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv('Stream_LR.csv',
                        header = TRUE) %>%
  select(From_Node, To_Node, Shape_Leng, UID, BARRIER_GROUP, Passability)

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
                                         EdgeType_field  = BARRIER_GROUP,
                                         Edge_attributes = c(UID, BARRIER_GROUP, Passability),
                                         Node_attributes = c(From_Node, Shape_Leng),
                                         Outlet_node     = outlet,
                                         graphFile       = 'Graph_LR.png')

# Extract graph object
graph.stream <- data.graph$data.graph

# Plot to confirm
gg0 <- ggnetwork(graph.stream,
                 layout =  layout_as_tree(graph.stream %>% as.undirected, root = 857),
                 scale = FALSE)

# windows()
plot <- 
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_nodes(alpha = 0.3) +
  geom_edges(alpha = 0.5,
             arrow = arrow(length = unit(10, "pt"), type = "closed"),
             aes(color = type)) +
  scale_color_viridis(discrete = TRUE)+
  geom_nodetext(aes(label = name), fontface = "bold") +
  theme_blank()
ggsave('Graph_TEST.png', plot = plot, 
       width = 30, height = 30, units = 'in',
       dpi = 800)


# Assign barrier passability  
# field.pass <- c('pass_u', 'pass_d')
# for(i in 1:length(field.pass)){
#   graph.stream <- set_edge_attr(graph.stream,
#                                 field.pass[i],
#                                 value = ifelse(E(graph.stream)$type == 'Confluence',
#                                                1.0,
#                                                0.5))
# }

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

# Calculate connectivty indices

# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity
out.dci <- index_calculation(graph = graph.stream,
                             weight = 'Shape_Leng',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------



View(as_data_frame(graph.stream, what = 'edges'))
View(as_data_frame(graph.stream, what = 'vertices'))
