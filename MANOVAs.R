# Title: MANOVAs

# Name: Heather White

# Date created: 01/09/21

# Last modified: 18/11/22

# License: MIT license


# MANOVA to explore correlates with skull shape 

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(ape)


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
# Access developmental strategy
Dev <- info$Dev_strategy


# Perform allometry and access residuals
fit1 <- procD.lm(absent_Proc~log(CS))
residuals <- fit1$residuals

# Create a geomorph dataframe
gdf <- geomorph.data.frame(absent_Proc, Age = info$Age, Species = info$Species, Dev_strat = info$Dev_strategy, Infraclass = info$Subclass, residuals = residuals)
# Name the dataframe
names(gdf) <-c("Proc_coords","Age","Species", "Development_strategy", "Infraclass","Residuals")


#######################################################################################################

# STEP 2: Perform MANOVAs on full dataset


# MANOVA for differences in skull shape based on age category
ageMVA <-procD.lm(Proc_coords ~ Age, f2 = NULL, f3 = NULL, logsz = TRUE, data = gdf, 
                  iter = 999, print.progress = FALSE)
summary(ageMVA)

# Residuals vs fitted plot
plot(ageMVA, method="PredLine")
# Looking at outliers
plot(ageMVA, outliers = TRUE) 

# Evaluate the assumptions of the model
par(mfrow = c(2,2))
plot(ageMVA)
par(mfrow = c(1,1))


######

# MANOVA for differences in skull shape based on developmental strategy
devMVA <-procD.lm(Residuals ~ Development_strategy, f2 = NULL, f3 = NULL, data = gdf, 
                  iter = 999, print.progress = FALSE)
summary(devMVA)


#######################################################################################################

# STEP 3: Perform phylogenetic MANOVA - adults only


# Load adult data

# Read in phylogeny
my_tree <- read.nexus("Data/my_mammal_tree.nexus")
# Load the adult only Procrusted data
load(file = "Data/adult_Procrusted_LMs.Rdata")
# Load adults only specimen details
adult_specimen_details <- read.csv("Data/Specimen_info_adults.csv", header = T, sep=",")
# Read in the specimen names to match the tree
tree_names <- read.csv("Data/tree_taxa_names.csv")
tree_names <- tree_names$Taxa_names

# Unname the adult LM dataset - required for later analysis
adults <- unname(adults)

# Check the tree tip names match the shape data
shape.names <-dimnames(adults)[[3]]
tree.names <-my_tree$tip.label
setdiff(shape.names,tree.names)
# Set the shape data to have the same names as the phylogeny
dimnames(adults)[[3]] <- tree_names


# Format the data

# Add the phylogeny to the geomorph dataframe
gdf_adults <- geomorph.data.frame(adults, phy = my_tree, Species = adult_specimen_details$Species,
                           Dev_strategy = adult_specimen_details$Precocial_altricial_spectrum,
                           CS = adult_specimen_details$CS)
# Name the geomorph dataframe 
names(gdf_adults) <-c("Phy","Proc_coords","Species", "Dev_strategy", "CS")


# Perform pMANOVA

# Phylogenetic MANOVA for developmental strategy
devPMVA <- procD.pgls(Proc_coords ~ Dev_strategy, phy = Phy, data = gdf_adults, iter = 999)                          
summary(devPMVA)


#######################################################################################################

# STEP 4: Perform MANOVA and pMANOVA for placentals only

# Subset out the placentals
placentals <- filter(info, Subclass == "Eutheria")
# Select the placental only LMs
plac_no <- placentals[,1]
placental_LMs <- absent_Proc[,,plac_no]

# Get the raw Procrustes data into a geomorph dataframe
gdf_plac <- geomorph.data.frame(placental_LMs, Age = placentals$Age, 
                                Species = placentals$Species, 
                                Dev_strat = placentals$Dev_strategy, 
                                Infraclass = placentals$Subclass)
# name first part of data frame (containing Procrustes data)
names(gdf_plac) <-c("Proc_coords","Age","Species", "Development_strategy", "Infraclass")

###

# Perform the MANOVA for differences in skull shape based on developmental strategy
devMVA_plac <-procD.lm(Proc_coords ~ Development_strategy, f2 = NULL, f3 = NULL, data = gdf_plac, 
                  iter = 999, print.progress = FALSE)
summary(devMVA_plac)

###

# Perform the pMANOVA

# Read in placental only phylogeny
my_tree_placentals <- read.nexus("Data/my_mammal_tree_placentals.nexus")

# Select the placentals only
adult_placental_LMs <- adults[,,-c(1,12,14,16,19,20,22)]
# Select the placental specimen details
adult_info_plac <- filter(adult_specimen_details, Subclass == "Eutheria")
# Read in the specimen names to match the tree
tree_names <- read.csv("Data/tree_taxa_names_placentals.csv")
tree_names <- tree_names$Taxa_names

# Unname the placental LMs for later analysis
adult_placental_LMs <- unname(adult_placental_LMs)
# Check the tree tip names match the shape data
shape.names <-dimnames(adult_placental_LMs)[[3]]
tree.names <-my_tree_placentals$tip.label
setdiff(shape.names,tree.names)
# Set the shape data to have the same names as the phylogeny
dimnames(adult_placental_LMs)[[3]] <- tree_names

# Create geomorph dataframe for placental adults
gdf_plac_adults <- geomorph.data.frame(adult_placental_LMs, phy = my_tree_placentals, 
                           Species = adult_info_plac$Species,
                           Dev_strategy = adult_info_plac$Precocial_altricial_spectrum)
# Name the geomorph dataframe 
names(gdf_plac_adults) <-c("Phy","Proc_coords","Species", "Dev_strategy")


# Perform pMANOVA for developmental strategy and skull shape
devPMVA_plac <- procD.pgls(Proc_coords ~ Dev_strategy, phy = Phy, data = gdf_plac_adults, iter = 999)                          
summary(devPMVA_plac)


### END ###
