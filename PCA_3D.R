# Title: 3D PCA

# Name: Heather White

# Date created: 31/08/21

# Last modified: 18/11/22

# License: MIT license


# Plot 3D PCA with developmental age as the z-axis, and PC1 and PC2 as the x- and y-axes

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(plotly)
library(viridis)
library(RColorBrewer)

#######################################################################################################

# STEP 1: Load the data

# Load in the LM data (mirrored, Procrusted, RHS deleted, missing LMs estimated and variably present LMs slid)
load(file ="Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")

# Read in the specimen details
specimen_details <- read.csv("Data/Specimen_info.csv", header = T, sep=",")

#######################################################################################################

# STEP 2: Perform PCA

# Perform PCA
PCA<-gm.prcomp(absent_Proc)
# Obtain PCA summary info
PCA_summary <- summary(PCA)
PCA_summary <- PCA_summary$PC.summary


# Associate specimen info with the PCA results table
# Convert PCA results into tibble format
PCAresults_all<-as_tibble(PCA$x) 
PCAresults_all<-PCAresults_all %>% mutate(Major_clade = specimen_details$Major_clades, 
                                            Species = specimen_details$Species, 
                                            Original_age = specimen_details$Age, 
                                            Percent_age = specimen_details$CS_percent_adult,
                                            Age_cat = specimen_details$CS_binned, 
                                            Percent_binned = specimen_details$Percent_binned)


#######################################################################################################

# STEP 3: Plot the PCAs


# Before plotting the PCAs - give the axes names
axx <- list(title = "PC1 (31.16%)")
axy <- list(title = "PC2 (14.12%)")
axz <- list(title = "Developmental age (%)")

# Plot PCA with PC1 and PC2 on x- and y-axes and developmental age (% of adult) on y-axis
P1 <- plot_ly(PCAresults_all, x=~Comp1, y=~Comp2, z=~Percent_age, type="scatter3d", 
              mode="markers", color = ~Percent_age)
P1 <- P1 %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz),
         title = "3D PCA with developmental age",
         legend=list(title=list(text='Developmental age (%)')))
P1


######

# Plot PCA with continuous and discrete age categories
P2 <- plot_ly(PCAresults_all, x=~Comp1, y=~Comp2, z=~Percent_age, type="scatter3d", 
              mode="markers", color = ~Percent_age, symbol = ~Age_cat, 
              symbols = c('circle', 'square-open', 'square', 'circle-open'), 
              marker = list(size =7))
P2 <- P2 %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz),
         title = "3D PCA with developmental age",
         legend=list(title=list(text='Developmental age (%)')))
P2

######

# Comparing discrete age category approaches
# Plot PCA with continuous age (% of adult) binned into 4 groups 
# Alongside the four discrete age categories (F, I, SA, A)


# Plot PCA with binned continuous age on z-axis and discrete age categories as colours/symbols
P3 <- plot_ly(PCAresults_all, x=~Comp1, y=~Comp2, z=~Percent_binned, 
                            type="scatter3d", mode="markers", 
                            color = ~Percent_binned, 
                            colors = c("#440154FF", "#33638DFF", "#3CBB75FF", "#FDE725FF"),
                            symbol = ~Age_cat, 
                            symbols = c('circle', 'square-open', 'square', 'circle-open'), 
                            marker = list(size =7))
P3 <- P3 %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz),
         title = "3D PCA with developmental age",
         legend=list(title=list(text='Developmental age (%)')),
         yaxis = list(categoryarray = ~Percent_binned, categoryorder = "category ascending"), scale_y_discrete(labels=c("Group1" = "Group 1 (7.37 - 30.6 %)", "Group2" = "Group 2 (30.6 - 53.7 %)", 
                                           "Group3" = "Group 3 = 53.7 - 76.9 %", "Group4" = "76.9 - 100 %")))
P3

### END ###
