# Title: Allometry

# Name: Heather White

# Date created: 25/08/21

# Last modified: 18/11/22

# License: MIT license


# Linear models looking at the relationship between Procrustes shape variables and CS - allometry

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)


#######################################################################################################

# STEP 1: Load and format the data

# Load the LM data
load(file = "Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")

# Load the specimen details
info <- read.csv(file = "Data/Specimen_info.csv", header =T, sep = ",")
# Access the specimen CS
CS <- info$CS
# Access the species
Species <- info$Species

# Remove specimen names from the coordinate dataset - need to be as numbers for later analysis
absent_Proc <- unname(absent_Proc)


#######################################################################################################

# STEP 2: Perform allometry

# Fit LM with CS - shape is correlated by size (to test allometry)
fit1 <- procD.lm(absent_Proc~log(CS)) # 1000 permutations
summary(fit1)

allometry_results <- fit1$aov.table

# Residuals vs fitted plot
plot(fit1, method="PredLine")
# Looking at outliers
plot(fit1, outliers = TRUE) 

# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(fit1)
par(mfrow = c(1,1))

#######

# Fit the model with CS and species
fit2 <- procD.lm(absent_Proc~log(CS)*Species) # 1000 permutations automatic; CS is logged
summary(fit2)

# Residuals vs fitted plot
plot(fit2, method="PredLine")
# Looking at outliers
plot(fit2, outliers = TRUE) # This line plots the outliers graph and highlights the outliers in red
# Here the Cebus and Ornithorhynchus are the outliers

# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(fit1)
par(mfrow = c(1,1))


##########################################################################################################

# STEP 3: Plot allometry results

# Produce geomorph dataframe with Proc_LMs
gdf <- geomorph.data.frame(absent_Proc, CSize = info$CS, Species = info$Species, Specimens = info$Specimen_name)
# Name the dataframe
names(gdf) <-c("Proc_coords","CSize","Species", "Specimens")

# Assigning species' colours
col.species = c("paleturquoise1", "pink1", "darkolivegreen1", "palevioletred",
                "lawngreen", "mediumvioletred", "burlywood1", "sandybrown", "mediumpurple1",
                "darkorange2", "green3", "turquoise1", "chartreuse4", "gold1",
                "orangered1", "cyan3", "darkgreen", "darkorchid4",
                "dodgerblue2", "blue2", "red3", "midnightblue")[gdf$Species]


# Plot allometry - coloured by species
plotAllometry(fit1, size = gdf$CSize, method = "CAC", pch = 16, col = col.species, cex = 1.5)
# CAC = common allometric component

# Check assumptions
plot(fit1, type = "regression", predictor = CS)
par(mfrow=c(2, 2), mar=c(5,4,4,2))
plot(fit1, type = "diagnostics")


### END ###
