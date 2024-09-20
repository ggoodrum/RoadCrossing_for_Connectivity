# Title: Connectivity_Testing
# Author: Greg Goodrum
# Last update: 7/14/2023
# Contact: greg.goodrum@usu.edu
# Description: Testing workflow for calculating stream network connectivity
#              Based on the following workflow:
#              https://damianobaldan.github.io/riverconn_tutorial/

# --------------------------------------------------------------------- #


# --------------------------------------------------------------------- #
# 01. Set up workspace
# ---------------------------------------------------------------------

# Summary: Set up workspace, load data, and load relevant packages.

# Clean workspace
rm(list = ls())

# Load Packages
if(!require("rstudioapi")){
  install.packages("rstudioapi", dependencies = T); library(rstudioapi)}
if(!require("dplyr")){
  install.packages("dplyr", dependencies = T); library(dplyr)}
if(!require("igraph")){
  install.packages("igraph", dependencies = T); library(igraph)}
if(!require("lubridate")){
  install.packages("lubridate", dependencies = T); library(lubridate)}
if(!require("foreach")){
  install.packages("foreach", dependencies = T); library(foreach)}
if(!require("parallel")){
  install.packages("parallel", dependencies = T); library(parallel)}
if(!require("snow")){
  install.packages("snow", dependencies = T); library(snow)}
if(!require("doSNOW")){
  install.packages("doSNOW", dependencies = T); library(doSNOW)}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 02. Initialize functions
# ---------------------------------------------------------------------

# FUNCTION: Prepare input network for DCI calculation
# https://github.com/cbedge/DCI/blob/master/networkprep_V2.R
calc.prep <- function(Parth, Quality = "Default"){
  Edge <- cbind(Parth$Fnode, Parth$Tnode, Parth$PermUS, Parth$PermDS) #The list of edges, change to now include US perm
  
  #Turning the list of edges into an adjacency matrix
  #For this to run correctly segments must be labelled 1, 2, 3... with no missing segments.
  Edge.adj <- matrix(0, nrow=length(Edge[,1]), ncol=length(Edge[,1]))
  for (i in 1:length(Edge[,1]))
  {
    Edge.adj[Edge[i,1], Edge[i,2]] <- Edge[i,3] #now a weighted matrix by US perm
    Edge.adj[Edge[i,2], Edge[i,1]] <- Edge[i,4] #weight by DS perm
  }
  
  #Now we create a directed graph of the river network from the adjacency matrix
  graph <- graph.adjacency(t(Edge.adj), mode="directed", weighted=TRUE) #Create a weighted graph
  
  #Now we create a Node file that contains all the attributes to run DCI
  Nodes <- data.frame(matrix(NA, nrow=length(Parth[,1]), ncol=6))
  colnames(Nodes) <- c("Junction", "ID", "Area", "Qual", "PermUS", "PermDS")
  
  Nodes$ID <- Parth$Fnode #upstream segment
  Nodes$Junction <- Parth$BARRIER_CO #barrier downstream of the segment, if present
  Nodes$Area <- Parth$Shape_Leng #length of segment
  
  if(Quality == "Default"){
    Nodes$Qual <- 1 #default value
  }else{
    Nodes$Qual <- Parth[,Quality]/100
  }
  Nodes$PermDS <- Parth$PermDS #permeability to move from focal segment to downstream segment
  Nodes$PermUS <- Parth$PermUS #assigning a permeability to move from segment downstream into the focal segment
  
  list("graph" = graph, "Nodes" = Nodes, "edgematrix" = Edge.adj)
}


# FUNCTION: Calculate DCIp and DCId
# https://github.com/cbedge/DCI/blob/master/DCI_V2.R
# NOTE: Edited for parallel by GG
DCI.calc <- function(nodefile, graphfile, mouth, edgemat, calc.DCIs = 0){
  #Construct a table with the start and end of each path between all segments, length of that path,
  
  # Set start time for progress tracking
  time.start <- as_datetime(Sys.time())
  
  # Delete previous log files
  file.remove('log.txt')
  
  # Detect cores
  parallel::detectCores()
  
  # Set number of cores, leaving one free for other tasks
  n.cores <- parallel::detectCores() - 1
  
  # Create cluster
  my.cluster <- snow::makeCluster(spec = n.cores,
                                  type = 'SOCK')
  
  # Register cluster
  doSNOW::registerDoSNOW(my.cluster)
  
  # Clear log before loop
  # writeLines(c(""), "log.txt")

  # Calculate c(i,j) parameter for all connections in network  
  path.vals <- foreach(k = 1:length(nodefile[,1]),
                       .combine = 'rbind',
                       .packages = c('igraph',
                                     'lubridate')) %dopar% {

                         # Calculate time at start of loop and model progression              
                         time.now  <- as_datetime(Sys.time())
                         time.tot  <- round(seconds_to_period(as.numeric(difftime(time.now, time.start, units = 'secs'))), digits = 2)
                         
                         # Write output to start of loop
                         log.text <- paste0("Starting node: ", k, " / ", length(nodefile[,1]), "  ",
                                            'Processing Time: ', time.tot)
                         write.table(log.text, "log.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
                         
                         # Determine all paths and calculate c(i,j)
                         path <- all_simple_paths(graph = graphfile, from = k, to = V(graphfile), mode = "all")
                         start.seg <- sapply(path, "[[", 1) #Identify the starting segment
                         end.seg <- sapply(path, tail, 1) #Identify the final segment
                         path.vals2 <- as.data.frame(cbind(start.seg, end.seg)) #combine start and end segments together
                         for (j in 1:length(path)) {
                           cij <- 1 #Default value for cij
                           for (i in 1:(length(path[[j]])-1)) {
                             cij <- cij * (edgemat[path[[j]][i], path[[j]][i+1]] * edgemat[path[[j]][i+1], path[[j]][i]])
                           }
                           path.vals2$cij[j] <- cij
                         }
                         itself <- c(k,k,  1)
                         path.vals <- rbind( path.vals2, itself)
                       }
  
  # Stop cluster
  snow::stopCluster(my.cluster)
  
  # Attach stream lengths to paths in network 
  path.vals$Start.Length <- (nodefile$Area[match(path.vals$start, nodefile$ID)]) * (nodefile$Qual[match(path.vals$start, nodefile$ID)])
  path.vals$End.Length <- (nodefile$Area[match(path.vals$end, nodefile$ID)]) * (nodefile$Qual[match(path.vals$end, nodefile$ID)])  
  
  #DCIp Calculation
  tot.length <- sum(nodefile$Area)
  DCIp <- 0
  for (i in 1:nrow(path.vals)) {
    DCIp <- DCIp + (path.vals$cij[i] * (path.vals$Start.Length[i]/tot.length) * (path.vals$End.Length[i]/tot.length))
  }  
  
  #DCId Calculation
  DCId.data <- path.vals[path.vals$end.seg==mouth, ]
  DCId.data$val1 <- DCId.data$cij * (DCId.data$Start.Length/tot.length)
  DCId <- sum(DCId.data$val1)
  
  #Combine DCIp and DCId into output
  DCI.results <- data.frame(DCIp, DCId)
  colnames(DCI.results) <- c("DCIp", "DCId")
  
  #Calculating DCIs, recalculate DCId for each segement as the mouth of the network
  if(calc.DCIs == 1){
    DCIs.vals <- data.frame(ID=integer(length(nodefile$ID)), DCIs=double(length(nodefile$ID)) )
    j <- 1
    for (i in unique(nodefile$ID)){
      DCId.data2 <- path.vals[path.vals$end.seg==i, ]
      DCId.data2$val1 <- DCId.data2$cij * (DCId.data2$Start.Length/tot.length)
      DCIs <- sum(DCId.data2$val1)
      DCIs.vals$ID[j] <- i
      DCIs.vals$DCIs[j] <- DCIs
      j <- j+1
    }
  }else{
    DCIs.vals <- "DCIs was not calculated. To caluclate DCIs set calc.DCIs to 1"
  }
  
  list("Index" = DCI.results, "DCIs" = DCIs.vals)
}

# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 03. Test Example
# ---------------------------------------------------------------------

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv('Stream_Test.csv',
                        header = TRUE)

# Add fields for DCI calculation
data.stream <- data.stream %>% mutate(Fnode      = From_Node,
                                      Tnode      = To_Node,
                                      BARRIER_CO = UID,
                                      PermUS     = 1,
                                      PermDS     = 1,
                                      Qual       = 100)

# Slim data data for testing
data.test <- data.stream %>% select(Fnode,
                                    Tnode,
                                    Shape_Leng,
                                    BARRIER_CO,
                                    PermUS,
                                    PermDS,
                                    Qual)

# # OPTION: Set all barriers to US/DS passability = 0.5
data.test <- data.test %>% mutate(PermUS = case_when(grepl('UID', BARRIER_CO) ~ 0.5,
                                                     .default = 1.0))

# Check that To_Nodes all come from From_Nodes
# NOTE: This should return 1 value, which should match to the outlet
data.outlet <- data.test %>% filter(!Tnode %in% Fnode)
# View(data.outlet)

# Check that nodes occur no more than twice in To_Node (binary junctions)
data.confluence <- as.data.frame(data.test %>% group_by(Tnode) %>% filter(n() > 2))

# Add terminal node
data.test <- data.test %>% add_row(Fnode      = data.outlet$Tnode,
                                   Tnode      = NA,
                                   Shape_Leng = 0,
                                   BARRIER_CO = NA,
                                   PermUS     = 1,
                                   PermDS     = 1,
                                   Qual       = 0)


# Prep data
data.prep <- calc.prep(Parth = data.test, Quality = "Default")

# Calculate DCI
data.dci <- DCI.calc(nodefile  = data.prep$Nodes, 
                     graphfile = data.prep$graph, 
                     edgemat   = data.prep$edgematrix, 
                     mouth     = data.outlet$Tnode, 
                     calc.DCIs = 0)
data.dci$Index
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 04. UBR Test
# ---------------------------------------------------------------------

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv('Stream_UBR.csv',
                        header = TRUE)

# Add fields for DCI calculation
data.stream <- data.stream %>% mutate(Fnode      = From_Node,
                                      Tnode      = To_Node,
                                      BARRIER_CO = UID,
                                      PermUS     = 1,
                                      PermDS     = 1,
                                      Qual       = 100)

# Slim data data for testing
data.test <- data.stream %>% select(Fnode,
                                    Tnode,
                                    Shape_Leng,
                                    BARRIER_CO,
                                    PermUS,
                                    PermDS,
                                    Qual)

# OPTION: Set all barriers to US/DS passability = 0.5
# data.test <- data.test %>% mutate(PermUS = case_when(grepl('UID', BARRIER_CO) ~ 0.5,
#                                                      .default = 1.0))

# Check that To_Nodes all come from From_Nodes
# NOTE: This should return 1 value, which should match to the outlet
data.outlet <- data.test %>% filter(!Tnode %in% Fnode)
# View(data.outlet)

# Check that nodes occur no more than twice in To_Node (binary junctions)
data.confluence <- as.data.frame(data.test %>% group_by(Tnode) %>% filter(n() > 2))
# View(data.confluence)

# Add terminal node
data.test <- data.test %>% add_row(Fnode      = data.outlet$Tnode,
                                   Tnode      = NA,
                                   Shape_Leng = 0,
                                   BARRIER_CO = NA,
                                   PermUS     = 1,
                                   PermDS     = 1,
                                   Qual       = 0)


# Prep data
data.prep <- calc.prep(Parth = data.test, Quality = "Default")

# Calculate DCI
data.dci <- DCI.calc(nodefile  = data.prep$Nodes, 
                     graphfile = data.prep$graph, 
                     edgemat   = data.prep$edgematrix, 
                     mouth     = data.outlet$Tnode, 
                     calc.DCIs = 0)
data.dci$Index
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 05. LR Test
# ---------------------------------------------------------------------

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv('Stream_LR.csv',
                        header = TRUE)

# Add fields for DCI calculation
data.stream <- data.stream %>% mutate(Fnode      = From_Node,
                                      Tnode      = To_Node,
                                      BARRIER_CO = UID,
                                      PermUS     = 1,
                                      PermDS     = 1,
                                      Qual       = 100)

# Slim data data for testing
data.test <- data.stream %>% select(Fnode,
                                    Tnode,
                                    Shape_Leng,
                                    BARRIER_CO,
                                    PermUS,
                                    PermDS,
                                    Qual)

# OPTION: Set all barriers to US/DS passability = 0.5
# data.test <- data.test %>% mutate(PermUS = case_when(grepl('UID', BARRIER_CO) ~ 0.5,
#                                                      .default = 1.0))

# Check that To_Nodes all come from From_Nodes
# NOTE: This should return 1 value, which should match to the outlet
data.outlet <- data.test %>% filter(!Tnode %in% Fnode)
# View(data.outlet)

# Check that nodes occur no more than twice in To_Node (binary junctions)
data.confluence <- as.data.frame(data.test %>% group_by(Tnode) %>% filter(n() > 2))
# View(data.confluence)

# Add terminal node
data.test <- data.test %>% add_row(Fnode      = data.outlet$Tnode,
                                   Tnode      = NA,
                                   Shape_Leng = 0,
                                   BARRIER_CO = NA,
                                   PermUS     = 1,
                                   PermDS     = 1,
                                   Qual       = 0)


# Prep data
data.prep <- calc.prep(Parth = data.test, Quality = "Default")

# Calculate DCI
data.dci <- DCI.calc(nodefile  = data.prep$Nodes, 
                     graphfile = data.prep$graph, 
                     edgemat   = data.prep$edgematrix, 
                     mouth     = data.outlet$Tnode, 
                     calc.DCIs = 0)
data.dci$Index
# ---------------------------------------------------------------------

# --------------------------------------------------------------------- #
# 06. BR Test
# ---------------------------------------------------------------------

# Declare working directory
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Load in data
data.stream <- read.csv('Stream_BR.csv',
                        header = TRUE)

# Add fields for DCI calculation
data.stream <- data.stream %>% mutate(Fnode      = From_Node,
                                      Tnode      = To_Node,
                                      BARRIER_CO = UID,
                                      PermUS     = 1,
                                      PermDS     = 1,
                                      Qual       = 100)

# Slim data data for testing
data.test <- data.stream %>% select(Fnode,
                                    Tnode,
                                    Shape_Leng,
                                    BARRIER_CO,
                                    PermUS,
                                    PermDS,
                                    Qual)

# OPTION: Set all barriers to US/DS passability = 0.5
# data.test <- data.test %>% mutate(PermUS = case_when(grepl('UID', BARRIER_CO) ~ 0.5,
#                                                      .default = 1.0))

# Check that To_Nodes all come from From_Nodes
# NOTE: This should return 1 value, which should match to the outlet
data.outlet <- data.test %>% filter(!Tnode %in% Fnode)
# View(data.outlet)

# Check that nodes occur no more than twice in To_Node (binary junctions)
data.confluence <- as.data.frame(data.test %>% group_by(Tnode) %>% filter(n() > 2))
# View(data.confluence)

# Add terminal node
data.test <- data.test %>% add_row(Fnode      = data.outlet$Tnode,
                                   Tnode      = NA,
                                   Shape_Leng = 0,
                                   BARRIER_CO = NA,
                                   PermUS     = 1,
                                   PermDS     = 1,
                                   Qual       = 0)


# Prep data
data.prep <- calc.prep(Parth = data.test, Quality = "Default")

# Calculate DCI
data.dci <- DCI.calc(nodefile  = data.prep$Nodes, 
                     graphfile = data.prep$graph, 
                     edgemat   = data.prep$edgematrix, 
                     mouth     = data.outlet$Tnode, 
                     calc.DCIs = 0)
data.dci$Index
# ---------------------------------------------------------------------

