# Title: RTK Data Processing
# Author: Greg Goodrum
# Last update: 17/06/2023
# Contact: greg.goodrum@usu.edu
# Description: Processing raw RTK data for input into ArcGIS


# --------------------------------------------------------------------- #
# 00. Instructions
# ---------------------------------------------------------------------

# Summary: This code reads in the raw ASCII RTK data and converts it to
#          a .csv for projection into a GIS.


# --------------------------------------------------------------------- #
# 01. Convert ASCII to CSV
# ---------------------------------------------------------------------

  # Summary: Set up workspace, load data, and load relevant packages.

  # Clean workspace
    rm(list = ls())

  # Load Packages
    if(!require("dplyr")){
      install.packages("dplyr", dependencies = T); library(dplyr)}
    if(!require("ggplot2")){
      install.packages("ggplot2", dependencies = T); library(ggplot2)}
    if(!require("rstudioapi")){
      install.packages("rstudioapi", dependencies = T); library(rstudioapi)}

  # Declare working directory
    pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
    setwd(pwd)
    
  # Load data
    data.2022 <- read.csv('BARRIER2022.txt', header = FALSE)
    data.2023 <- read.csv('BARRIER2023.txt', header = FALSE)
    
  # Add sample year column
    data.2022 <- data.2022 %>% mutate(Year = '2022')
    data.2023 <- data.2023 %>% mutate(Year = '2023')
    
  # Combine datasets
    data.RTK <- rbind(data.2022, data.2023)
    
  # Remove and rename columns
    data.RTK <- data.RTK %>% select(V1, V2, V3, V4, Year) %>%
                             rename(Site_ID = 'V1',
                                    Y_Northing_UTM = 'V2',
                                    X_Easting_UTM = 'V3',
                                    Elevation_m = 'V4')
    
  # Export data
    write.csv(data.RTK,
              file = 'Barrier_RTK_2022_2023.csv', row.names = FALSE)
    
# ---------------------------------------------------------------------