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
# 02. Test Example
# ---------------------------------------------------------------------

# Clean workspace
rm(list = ls())

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Pre-process data ----------------------------------------------------

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv('Stream_Test.csv',
                        header = TRUE) %>%
               select(From_Node, To_Node, Shape_Leng, UID, BARRIER_GROUP, Passability)

# Check that From_Nodes are all unique
check.reaches    <- as.data.frame(data.stream %>% group_by(From_Node) %>% filter(n() > 1))
nrow(check.reaches) # Check if value != 0

# Check that nodes occur no more than twice in To_Node (binary junctions)
check.confluence <- as.data.frame(data.stream %>% group_by(To_Node) %>% filter(n() > 2))
nrow(check.confluence) # Check if value != 0

# Check that To_Nodes all come from From_Nodes
# NOTE: This should return 1 value, which should match to the outlet
data.outlet <- data.stream %>% filter(!To_Node %in% From_Node)
outlet <- data.outlet$From_Node # Check if more than 1 record


# Create igraph object ------------------------------------------------

# Remove row with terminal node and edge attributes
data.graph <- data.stream %>% select(From_Node, To_Node, UID, BARRIER_GROUP, Passability) %>%
                              filter(To_Node != data.outlet$To_Node)

# Create igraph object with edge attributes
graph.stream <- igraph::graph_from_data_frame(data.graph)

# Attribute vertices with stream reach data
data.vertices <- data.stream %>% select(From_Node, Shape_Leng)
graph.stream  <- set_vertex_attr(graph.stream, 'length_km', value = data.vertices$Shape_Leng)

# Check vertex attribution
check.vertices <- igraph::as_data_frame(graph.stream, what = 'vertices') %>%
                  mutate(name = as.numeric(name))
check.verts    <- left_join(check.vertices, data.vertices, by = c('name' = 'From_Node'))
check.vdupl    <- check.verts %>% group_by_all() %>%
                                  filter(n() >  1) %>%
                                  ungroup()
nrow(check.vdupl)

# Assign river network directionality
graph.stream <- set_graph_directionality(graph.stream, 
                                         field_name = 'name',
                                         outlet_name = outlet)

# Plot to confirm
# gg0 <- ggnetwork(graph.stream,
#                  layout =  layout_as_tree(graph.stream %>% as.undirected,
#                                           root = outlet), scale = FALSE)
# 
# windows()
# ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_nodes(alpha = 0.3) +
#   geom_edges(alpha = 0.5,
#              arrow = arrow(length = unit(10, "pt"), type = "closed")) +
#   geom_nodetext(aes(label = name), fontface = "bold") +
#   theme_blank()

# Add barrier passability fields ----------------------------------

# NOTE: Barriers are an edge attribute and represent either barriers or confluences

# Assign edge type as either a confluence OR barrier
graph.stream <- set_edge_attr(graph.stream,
                              'type',
                              value = ifelse(E(graph.stream)$BARRIER_GROUP == "",
                                             'Confluence', E(graph.stream)$BARRIER_GROUP))

# Initialize upstream and downstream passability fields
field.pass <- c('pass_u', 'pass_d')
for(i in 1:length(field.pass)){
  graph.stream <- set_edge_attr(graph.stream,
                                field.pass[i],
                                value = 1.0)
}

# Plot to confirm
# gg0 <- ggnetwork(graph.stream, 
#                  layout =  layout_as_tree(graph.stream %>% as.undirected, root = outlet), 
#                  scale = FALSE)
# 
# windows()
# ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_nodes(alpha = 0.3) +
#   geom_edges(alpha = 0.5, 
#              arrow = arrow(length = unit(10, "pt"), type = "closed"), 
#              aes(color = type)) + 
#   scale_color_viridis(discrete = TRUE)+
#   geom_nodetext(aes(label = name), fontface = "bold") +
#   theme_blank()

# Assign barrier passability  -----------------------------------------

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

# Calculate connectivty indices  --------------------------------------

# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity
out.dci <- index_calculation(graph = graph.stream,
                             weight = 'length_km',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 04. UBR Example
# ---------------------------------------------------------------------

# Clean workspace
rm(list = ls())

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Pre-process data ----------------------------------------------------

# Load in data
# NOTE: Input data represents the edge list
data.stream <- read.csv('Stream_UBR.csv',
                        header = TRUE) %>%
  select(From_Node, To_Node, Shape_Leng, UID, BARRIER_GROUP, Passability)

# Check that From_Nodes are all unique
check.reaches    <- as.data.frame(data.stream %>% group_by(From_Node) %>% filter(n() > 1))
nrow(check.reaches) # Check if value != 0

# Check that nodes occur no more than twice in To_Node (binary junctions)
check.confluence <- as.data.frame(data.stream %>% group_by(To_Node) %>% filter(n() > 2))
nrow(check.confluence) # Check if value != 0

# Check that To_Nodes all come from From_Nodes
# NOTE: This should return 1 value, which should match to the outlet
data.outlet <- data.stream %>% filter(!To_Node %in% From_Node)
outlet <- data.outlet$From_Node # Check if more than 1 record


# Create igraph object ------------------------------------------------

# Remove row with terminal node and edge attributes
data.graph <- data.stream %>% select(From_Node, To_Node, UID, BARRIER_GROUP, Passability) %>%
  filter(To_Node != data.outlet$To_Node)

# Create igraph object with edge attributes
graph.stream <- igraph::graph_from_data_frame(data.graph)

# Attribute vertices with stream reach data
data.vertices <- data.stream %>% select(From_Node, Shape_Leng)
graph.stream  <- set_vertex_attr(graph.stream, 'length_km', value = data.vertices$Shape_Leng)

# Check vertex attribution
check.vertices <- igraph::as_data_frame(graph.stream, what = 'vertices') %>%
  mutate(name = as.numeric(name))
check.verts    <- left_join(check.vertices, data.vertices, by = c('name' = 'From_Node'))
check.vdupl    <- check.verts %>% group_by_all() %>%
  filter(n() >  1) %>%
  ungroup()
nrow(check.vdupl)

# Assign river network directionality
graph.stream <- set_graph_directionality(graph.stream, 
                                         field_name = 'name',
                                         outlet_name = outlet)

# Plot to confirm
# gg0 <- ggnetwork(graph.stream,
#                  layout =  layout_as_tree(graph.stream %>% as.undirected,
#                                           root = outlet), scale = FALSE)
# 
# windows()
# ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_nodes(alpha = 0.3) +
#   geom_edges(alpha = 0.5,
#              arrow = arrow(length = unit(10, "pt"), type = "closed")) +
#   geom_nodetext(aes(label = name), fontface = "bold") +
#   theme_blank()

# Add barrier passability fields ----------------------------------

# NOTE: Barriers are an edge attribute and represent either barriers or confluences

# Assign edge type as either a confluence OR barrier
graph.stream <- set_edge_attr(graph.stream,
                              'type',
                              value = ifelse(E(graph.stream)$BARRIER_GROUP == "",
                                             'Confluence', E(graph.stream)$BARRIER_GROUP))

# Initialize upstream and downstream passability fields
field.pass <- c('pass_u', 'pass_d')
for(i in 1:length(field.pass)){
  graph.stream <- set_edge_attr(graph.stream,
                                field.pass[i],
                                value = 1.0)
}

# Plot to confirm
# gg0 <- ggnetwork(graph.stream, 
#                  layout =  layout_as_tree(graph.stream %>% as.undirected, root = outlet), 
#                  scale = FALSE)
# 
# windows()
# ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_nodes(alpha = 0.3) +
#   geom_edges(alpha = 0.5, 
#              arrow = arrow(length = unit(10, "pt"), type = "closed"), 
#              aes(color = type)) + 
#   scale_color_viridis(discrete = TRUE)+
#   geom_nodetext(aes(label = name), fontface = "bold") +
#   theme_blank()

# Assign barrier passability  -----------------------------------------

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

# Calculate connectivty indices  --------------------------------------

# NOTE: See Baldan et al. (2022) 'Introducing 'riverconn': an R package to assess river
#       connectivity indices' for information describing indices and calculations
#       https://www.sciencedirect.com/science/article/pii/S1364815222001748

# Calculate DCIp connectivity
out.dci <- index_calculation(graph = graph.stream,
                             weight = 'length_km',
                             B_ij_flag = FALSE,
                             dir_fragmentation_type = 'symmetric')

# ---------------------------------------------------------------------



View(as_data_frame(graph.stream, what = 'edges'))
View(as_data_frame(graph.stream, what = 'vertices'))
