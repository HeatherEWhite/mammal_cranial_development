# Title: Phylomorphospace

# Name: Heather White

# Date created: 07/10/21

# Last modified: 18/11/22

# License: MIT license


# Phylomorphospace - plotting adult mammal specimens only

#######################################################################################################

rm(list = ls())

library(ggplot2)
library(ggConvexHull)
library(ggpubr)
library(RColorBrewer)
library(ape)
library(geomorph)
library(tidyverse)
library(gginnards)
library(ggphylomorpho)

#######################################################################################################

# STEP 1: Load the data

# Read in the trimmed phylogeny
my_tree <- read.nexus("Data/my_mammal_tree.nexus")
# Read in taxa names that match the phylogeny
tree_names <- read.csv("Data/tree_taxa_names.csv")
tree_names <- tree_names$Taxa_names

# Load in the LM data (mirrored, Procrusted, RHS deleted, missing LMs estimated and variably present LMs slid)
load(file ="Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")

# Read in the specimen details that match the dataset
specimen_details <- read.csv("Data/Specimen_info.csv", header = T, sep=",")
adult_specimen_details <- read.csv("Data/Specimen_info_adults.csv", header = T, sep=",")


#######################################################################################################

# STEP 2: Plot the phylogeny

# Check the phylogeny with timescale 
plot(my_tree, cex = 0.5) # plotted in normal linear tree fashion
axisPhylo()


#######################################################################################################

# STEP 3: Format the data

# Seperate out the adults only from the Procrusted LM data into a new dataset
adults_Proc_coords <- absent_Proc[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]

# Check the tree tip names match the shape data
shape.names <-dimnames(adults_Proc_coords)[[3]]
tree.names <-my_tree$tip.label
setdiff(shape.names,tree.names)
# Set the shape data to have the same names as the phylogeny
dimnames(adults_Proc_coords)[[3]] <- tree_names


#######################################################################################################

# STEP 4: Perform PCA

# Perform PCA
PCA<-gm.prcomp(adults_Proc_coords)
# Obtain the proportion of variance for each PC
PCA_summary <- summary(PCA) 
PCA_summary <- PCA_summary$PC.summary
# Obtain the PC scores for each species for each PC
PCscores <- PCA$x

# Convert PCA results to a tibble
# Associate specimen info with the PCA
PCAresults<-as_tibble(PCA$x) # Converting the PCA results into tibble format, so we can access results for each PC to plot PCA in ggplot
PCAresults<-PCAresults %>% mutate(Species = adult_specimen_details$Species, Phylogeny_names = tree_names, Subclass = adult_specimen_details$Subclass, Major_clade = adult_specimen_details$Major_clades) # Add columns to the above table using the information from the spec_info file. This will help with grouping things when plotting the PCA

# Add tree name information to the tibble
pcscores_adults <- PCAresults %>%  mutate(.,ID = tree_names)
# Sort the table by custom Clade order 
pcscores_adults$Major_clade <- factor(pcscores_adults$Major_clade, levels = c("Afrotheria","Euarchontoglires", "Laurasiatheria", "Marsupialia", "Monotremata","Xenarthra"))


#######################################################################################################

# STEP 5: Plot the phylomorphospace

# Create a palette for the clades
my_palette_clade <- c("mediumpurple3", # Afrotheria
                      "#A6D854", # Euarchontoglires
                      "sandybrown", # Laurasiatheria
                      "cornflowerblue", # Marsupialia
                      "#FFD92F", # Monotremata
                      "palevioletred") # Xenarthra
image(1:6, 1, as.matrix(1:6), col = my_palette_clade, xlab = "Clade",
      ylab = "", yaxt = "n")



# Plot the tree and tip data
g <- ggphylomorpho(tree = my_tree, tipinfo = pcscores_adults, xvar=Comp1, yvar = Comp2, factorvar = Major_clade,labelvar = ID, tree.alpha = 0.7, edge.width = 0.5)
# Remove specimen IDs and point data - to make pretty in next step
g <- g %>%
  delete_layers("GeomPoint") %>%
  delete_layers("GeomTextRepel")
# Add information to make the plot pretty
g + geom_convexhull(data = PCAresults, aes(x=Comp1, y = Comp2, colour = Major_clade, fill = Major_clade), alpha = 0.4, show.legend = F) +
  geom_point(data = PCAresults, aes(x=Comp1, y = Comp2,colour = Major_clade, fill = Major_clade, shape = Subclass),size=3)+#,show.legend = FALSE)
  geom_text_repel(data = pcscores_adults, aes(x=Comp1, y=Comp2, label=Species, hjust=0,vjust=0, fontface ="italic"), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = my_palette_clade) +
  scale_fill_manual(values = my_palette_clade) +
  ggtitle("Phylomorphospace - PC1 and PC2") +
  guides(color=guide_legend("Clade"), fill = FALSE) +
  xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")) + 
  ylab(paste0("PC2 (",signif((PCA_summary[2,2]*100),4),"% of total variance)")) +
  xlim(-0.2,0.2) +
  theme_classic(base_size = 12)


### END ###
