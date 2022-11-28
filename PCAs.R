# Title: Principal Components Analysis

# Name: Heather White

# Date created: 20/08/21

# Last modified: 18/11/22

# License: MIT license


# Principal components analysis of the mirrored, Procrusted, and missing/absent LMs dealt with - LHS only

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(Morpho)
library(rgl)
library(ape)
library(paleomorph)
library(RRPP)
library(arrayhelpers)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(devtools)
library(ggConvexHull)
library(ggthemes)
library(ggfortify)
library(ggphylomorpho)
library(gginnards)
library(wesanderson)


#######################################################################################################

# STEP 1: Load and format the data

# Load the data that has been mirrored, Procrusted, RHS deleted, missing LMs estimated and variably present LMs slid
load(file ="Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")

# Load the specimen details
# Full dataset details:
specimen_details <- read.csv("Data/Specimen_info.csv", header = T, sep=",")
# Adults details only:
adult_specimen_details <- read.csv("Data/Specimen_info_adults.csv", header = T, sep=",")


# Seperate out groups of specimens from the LM data
# Adults only
adults <- absent_Proc[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]
save(adults, file = "Data/adult_Procrusted_LMs.Rdata")


####################################################################################################

# STEP 2: Perform PCA for the adults only

# Perform PCA
PCA_adults<-gm.prcomp(adults)
# Get a summary of all the PCA and PC axes
PCA_adults_summary <- summary(PCA_adults)
PCA_adults_summary <- PCA_adults_summary$PC.summary
# Access PC scores for each specimen
PCscores_adults <- PCA_adults$x

# Quickly plot the PCA
PCA.plot <-plot(PCA_adults, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA_adults$x[,1], y=PCA_adults$x[,2], rownames(PCA_adults$x))


# See 'PC_extreme_shapes.R' code for determining the shape variation for each PC

# Quickly view the shape variation for PC1

# Calculate the mean shape
shape <- mshape(adults)
# PC1 min shape
plotRefToTarget(shape, PCA_adults$shapes$shapes.comp1$min)
# PC1 max shape
plotRefToTarget(shape, PCA_adults$shapes$shapes.comp1$max)
# PC1 min shape
spheres3d(PCA_adults$shapes$shapes.comp1$min, radius = 0.002)
# PC1 max shape
spheres3d(PCA_adults$shapes$shapes.comp1$max, radius = 0.002)


# Adjust PCA data into format necessary for plotting

# Converting the PCA results into tibble format
PCAresults_adults<-as_tibble(PCA_adults$x) 
# Associate specimen info with the PCA
PCAresults_adults<-PCAresults_adults %>% mutate(Species = adult_specimen_details$Species, 
                                                Subclass = adult_specimen_details$Subclass, 
                                                Major_clade = adult_specimen_details$Major_clades, 
                                                Dev_strategy = adult_specimen_details$Precocial_altricial_spectrum) 


####################################################################################################

# STEP 3: Perform PCA for the full developmental dataset

# Perform PCA
PCA_all<-gm.prcomp(absent_Proc)
# Get a summary of all the PCA and PC axes
PCA_all_summary <- summary(PCA_all)
PCA_all_summary <- PCA_all_summary$PC.summary
# Access PC scores for each specimen
PCscores_all <- PCA_all$x

# Quickly plot the PCA
PCA.plot <-plot(PCA_all, cex = 2, pch = 16, axis1 = 1, axis2 = 2)
# To add specimen names to PCA
text(x=PCA_all$x[,1], y=PCA_all$x[,2], rownames(PCA_all$x))


# See 'PC_extreme_shapes.R' code for determining the shape variation for each PC

# Quickly view the shape variation for PC1

# Calculate the mean shape
shape <- mshape(absent_Proc)
# PC1 min shape
plotRefToTarget(shape, PCA_all$shapes$shapes.comp1$min)
# PC1 max shape
plotRefToTarget(shape, PCA_all$shapes$shapes.comp1$max)
# PC1 min shape
spheres3d(PCA_all$shapes$shapes.comp1$min, radius = 0.002)
# PC1 max shape
spheres3d(PCA_all$shapes$shapes.comp1$max, radius = 0.002)


# Adjust PCA data into format necessary for plotting

# Converting the PCA results into tibble format
PCAresults_all<-as_tibble(PCA_all$x) 
# Associate specimen info with the PCA
PCAresults_all<-PCAresults_all %>% mutate(Subclass = specimen_details$Subclass, 
                                          Major_clade = specimen_details$Major_clades, 
                                          Species = specimen_details$Species, 
                                          Dev_strategy = specimen_details$Dev_strategy, 
                                          Age = specimen_details$CS_binned_group, 
                                          Specimen = specimen_details$Specimen_name,
                                          CS_logged = specimen_details$CS_logged,
                                          CS_percent_adult = specimen_details$CS_percent_adult) 


# Multiply PC1 by -1 to flip the PC1 axis
# This is done to have youngest specimens falling at the negative of PC1 - more intuitive
PCAresults_PC1 <- PCAresults_all$Comp1 * -1

# Combine the rest of the PCA results the flipped PC1 results
PCAresults_all <- select(PCAresults_all, 2:172)
PCAresults_all <- cbind(PCAresults_PC1, PCAresults_all)
names(PCAresults_all)[names(PCAresults_all) == 'PCAresults_PC1'] <- 'Comp1'

# Save PC scores with specimen info as an R object
save(PC_scores_PC1_flipped, file = "Data/PC_scores_all_specimens.Rdata")

# Can view PC1 shapes again after this manipulation


# Order the data by clade and convert to factor for plotting and plot legend
PCAresults_all$Major_clade <- as.factor(PCAresults_all$Major_clade)
PCAresults_all <- PCAresults_all %>% arrange(Major_clade)
PCAresults_all

# Order the species to match the colour order below
PCAresults_all$Species <- factor(PCAresults_all$Species,
                                 levels = c("Macroscelides proboscideus", "Setifer setosus", "Cebus apella",            
                                            "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                                            "Rattus rattus", "Epomops franqueti", "Felis catus",       
                                            "Manis tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                                            "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                                            "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                                            "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                                            "Dasypus novemcinctus"),
                                 ordered = TRUE)


####################################################################################################

# STEP 4: Setup the colours for the PCA


# Create a colour palette fir the species
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                       "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                       "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                       "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:22, 1, as.matrix(1:22), col = mypalette_species, xlab = "Species",
      ylab = "", yaxt = "n")


# Colours for each major clade
col.clade = c("mediumpurple3", # Afrotheria
              "#A6D854", # Euarchontoglires
              "sandybrown", # Laurasiatheria
              "cornflowerblue", # Marsupialia
              "#FFD92F", # Monotremata
              "palevioletred") # Xenarthra
image(1:6, 1, as.matrix(1:6), col = col.clade, xlab = "Clade",
      ylab = "", yaxt = "n")


# Colours for developmental strategy
col.dev <- c("sandybrown", "dodgerblue4", "tan4", "lightskyblue3", "bisque")
image(1:5, 1, as.matrix(1:5), col = col.dev, xlab = "Dev_strategy",
      ylab = "", yaxt = "n")

# Colours for developmental age
col.age <- c("#440154FF", "#33638DFF", "#3CBB75FF", "#FDE725FF")
image(1:4, 1, as.matrix(1:4), col = col.age, xlab = "Age",
      ylab = "", yaxt = "n")


#######################################################################################################################################

# STEP 5: Plot PCAs


# Adult PCA - species labelled
PCA_adults_labelled <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, label=Species))+ 
  geom_point(size=1)+
  geom_text_repel(aes(fontface="italic"), size = 4)+
  theme_classic(base_size = 12)+ 
  xlab(paste0("PC3 (", signif((PCA_adults_summary[2,3]*100),4), "% of total variance)")) + 
  ylab(paste0("PC4 (",signif((PCA_adults_summary[2,4]*100),4),"% of total variance)")) +
  ggtitle("PCA adults labelled")
PCA_adults_labelled


# Adult PCA - developmental strategy
PCA_adults_dev_strategy <- ggplot(PCAresults_adults, aes(x=Comp1, y=Comp2, colour=Dev_strategy, fill = Dev_strategy))+ 
  geom_convexhull(alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(data = PCAresults_adults, aes(x=Comp1, y=Comp2, shape = Subclass , colour = Dev_strategy), size=3)+
  theme_classic(base_size = 12)+ 
  xlab(paste0("PC1 (", signif((PCA_adults_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_adults_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = col.dev)+
  scale_fill_manual(values = col.dev)+
  labs(color = "Developmental Strategy")+
  ggtitle("PCA adults divided by developmental strategy")
PCA_adults_dev_strategy


# Full dataset PCA - species labelled
PCA_all <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = Species, shape = Subclass))+ 
  geom_point(size=3)+
  theme_classic(base_size = 14)+ 
  xlab(paste0("PC1 (", signif((PCA_all_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_all_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = mypalette_species)+
  labs(color = "Clade")+
  ggtitle("PCA all species coloured")
PCA_all


# Full dataset PCA - discrete developmental age
PCA_all_age <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour=Age, fill = Age))+ 
  geom_convexhull(data = PCAresults_all, aes(x=Comp1, y = Comp2, fill = Age), alpha = 0.4, show.legend = FALSE, size =0.2)+
  geom_point(data = PCAresults_all, aes(x=Comp1, y=Comp2, shape = Subclass , colour = Age), size=2.5)+
  theme_classic(base_size = 14)+ 
  xlab(paste0("PC1 (", signif((PCA_all_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_all_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_manual(values = col.age, labels=c("Group1" = "Fetal", "Group2" = "Infant", 
                                                "Group3" = "Sub-adult", "Group4" = "Adult"))+
  scale_fill_manual(values = col.age, labels=c("Group1" = "Fetal", "Group2" = "Infant", 
                                               "Group3" = "Sub-adult", "Group4" = "Adult"))+
  labs(color = "Age")+
  ggtitle("PCA all - divided by age category")
PCA_all_age


#########

# Comparison between the two continuous developmental age metrics

# Full dataset PCA, where age = logged CS
PCA_all_CS <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = CS_logged))+ 
  geom_point(size=3)+
  theme_classic(base_size = 14)+ 
  xlab(paste0("PC1 (", signif((PCA_all_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_all_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_viridis()+
  ggtitle("PCA all - CS across all specimens")
PCA_all_CS


# Full dataset PCA, where age = percentage of adult CS 
PCA_all_CS_percent_adult <- ggplot(PCAresults_all, aes(x=Comp1, y=Comp2, colour = CS_percent_adult))+ 
  geom_point(size=3)+
  theme_classic(base_size = 14)+ 
  xlab(paste0("PC1 (", signif((PCA_all_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_all_summary[2,2]*100),4),"% of total variance)")) +
  scale_color_viridis()+
  ggtitle("PCA all - CS (percent adult) across all specimens")
PCA_all_CS_percent_adult


### END ###
