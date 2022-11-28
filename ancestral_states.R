# Title: Ancestral states allometry

# Name: Heather White

# Date created: 01/10/21

# Last modified: 18/11/22

# License: MIT license


# Ancestral state reconstruction from species ontogenetic allometric trajectory analysis
# Adapted from Morris et al. 2019

#######################################################################################################

rm(list = ls())

library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)
library(ggfortify)
library(RColorBrewer) 
library(borealis)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(png)
library(gridExtra)
library(phytools)
library(evomap)
library(abind)
library(dplyr)

require(devtools)
devtools::install_github("JeroenSmaers/evomap")
devtools::install_github("wabarr/ggphylomorpho")
devtools::install_github("aphanotus/borealis")

####################################################################################################

# STEP 1: Load the data

# Import trees in Nexus format - branch lengths needed
tree <- "Data/my_mammal_tree.nexus"
# Read the tree for analysis
tree_species <- read.nexus(tree) 
plot(tree_species)
# Check names 
summary(tree_species)

# Load the PC score data, procrusted coordinate data, and specimen information
load(file ="Data/PCAresults_all.Rdata")
load(file = "Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")
info <- read.csv(file = "Data/Specimen_info.csv", header =T, sep = ",")
Phylogeny_species_list <- read.csv("Data/tree_taxa_names.csv")


#######################################################################################################

# STEP 2: Format the data

# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
absent_Proc <- unname(absent_Proc)

# Combine data into geomorph dataframe
gdf <- geomorph.data.frame(absent_Proc, CSize = info$CS, logCSize = info$CS_logged, Species = info$Species, Specimens = info$Specimen_name)
# Name first part of data frame (containing Procrustes data)
names(gdf) <-c("Proc_coords","CSize", "logCSize", "Species", "Specimens")

# Pull out the species names info
species_list <- levels(as.factor(gdf$Species))


# Seperate the dataset by species
# Create a list containing a list for each species in the dataset
species_rows <- list()
for (i in 1:length(species_list)){
  species_rows[[i]] <- which(PCAresults_all$Species == species_list[i])
}

# Seperate out the shape (coords) and size data by species
# A list for each species containing this information is produced
coords_list <- list()
logCsize_list <- list()
for (i in 1:length(species_list)){
  coords_list[[i]] <- gdf$Proc_coords[,,species_rows[[i]]]
  logCsize_list[[i]] <- gdf$logCSize[species_rows[[i]]]
}

# Assocaite species names with the subsetted dataset above
# This has to be the same names as the phylogeny - must match
names(coords_list) = Phylogeny_species_list$Taxa_names
names(logCsize_list) = Phylogeny_species_list$Taxa_names


#####################################################################################################

# STEP 3: Perform allometry

# Calculate allometric regression for each species seperately 
allometry_species <- list()
for (i in 1:22){
  allometry_species[[i]] <- procD.lm(coords_list[[i]] ~ logCsize_list[[i]], iter = 999)
}

# View the allometry results - seperate into a list for each species
summary_allometry_species <- list()
for (i in 1:22){
  summary_allometry_species[[i]] <- summary(allometry_species[[i]])
 }

# Assign species names to the allometry results
# These have to match the phylogeny
names(summary_allometry_species) <- Phylogeny_species_list$Taxa_names
summary_allometry_species


######
# These results above are the same as the results from the allometry.R code
######


# Plot allometry to obtain regscores for allometry
allometry_species_plot <- list()
for (i in 1:22){
  allometry_species_plot[[i]] <- plot(allometry_species[[i]], type = "regression", predictor = logCsize_list[[i]], reg.type = "RegScore",
                                      main = "Shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)
}

# Pool the regression scores (shape) (y) and logCS (x) from the allometry plots above into a dataframe 
# Calculated for each species individually
x <- list()
for (i in 1:22){
x[[i]] <- c(allometry_species_plot[[i]][["plot.args"]][["y"]])
allometry_regscores <- unlist(x)
allometry_plot_data <- cbind(logCsize = info$CS_logged, RegScores = allometry_regscores)
}

# Convert dataframe to tibble and add necessary information
allometry_plot_data <- as_tibble(allometry_plot_data)
allometry_plot_data <- allometry_plot_data %>% mutate(Specimen_names = info$Specimen_name, Species = info$Species, 
                            Clade = info$Major_clades)
glimpse(allometry_plot_data)


###################################################################################################################

# STEP 4: Perform linear model to get slope/intercept for each species for allometry analysis

# Linear model for regression score (y) against logCS (x) 
# To get the line of best fit - slope and intercept
allometry_species_regline <- list()
for (i in 1:22){
  allometry_species_regline[[i]] <- lm(allometry_species_plot[[i]][["plot.args"]][["y"]] ~ allometry_species_plot[[i]][["plot.args"]][["x"]])
}

# Save the slope and intercept values from the linear model
SlopesList <- matrix()
InterceptsList <- matrix()
for (i in 1:22){
  SlopesList[[i]] <- allometry_species_regline[[i]][["coefficients"]][2]
  InterceptsList[[i]] <- allometry_species_regline[[i]][["coefficients"]][1]
}

# Assign species names to the allometry results
# These have to match the phylogeny
names(SlopesList) <- Phylogeny_species_list$Taxa_names
names(InterceptsList) <- Phylogeny_species_list$Taxa_names


########################################################################################################

# STEP 5: Calculate slopes and intercepts for ancestral states

# Calculate ancestral states using phytools for slope and intercept values
slope_anc.ML <- anc.ML(tree_species, SlopesList, CI = F)
int_anc.ML <- anc.ML(tree_species, InterceptsList, CI = F)

# Label the nodes of the phylogeny
plot(tree_species, cex=0.5, show.node.label = T)
nodelabels(cex=0.5)
# Save the ancestral node numbers
nodes <- length(slope_anc.ML$ace)
# Create col names for ancestral data frame
anc_col_names <- c("Slope","Intercept")

# Create dataframe with ancestral node values only
# Dataframe contains both slope and intercept values
anc_values_traj <- matrix(nrow = nodes, ncol = 2, dimnames = list(names(slope_anc.ML$ace),anc_col_names))
anc_values_traj[,"Slope"] <- slope_anc.ML$ace
anc_values_traj[,"Intercept"] <- int_anc.ML$ace


#######################################################################################################

# STEP 6: Format data for plotting

# Create dataframes for slope and intercept values
# Include both extant species values and calculated ancestral states
slope <- data.frame(Slope = c(SlopesList,anc_values_traj[,1]))
intercept <- data.frame(Intercept = c(InterceptsList,anc_values_traj[,2]))
# Combine the slope and intercept data into one dataframe for plotting
values_traj <- cbind(intercept,slope)

# Name the relevant ancestral clade nodes
rownames(values_traj)[rownames(values_traj) == "23"] <- "Ancestral mammal"
rownames(values_traj)[rownames(values_traj) == "24"] <- "Ancestral therian mammal"
rownames(values_traj)[rownames(values_traj) == "25"] <- "Ancestral placental mammal"
rownames(values_traj)[rownames(values_traj) == "27"] <- "Ancestral Laurasiatheria"
rownames(values_traj)[rownames(values_traj) == "31"] <- "Ancestral Euarchontoglires"
rownames(values_traj)[rownames(values_traj) == "36"] <- "Ancestral Afrotheria"
rownames(values_traj)[rownames(values_traj) == "37"] <- "Ancestral Xenarthra"
rownames(values_traj)[rownames(values_traj) == "39"] <- "Ancestral Marsupialia"

# Make tibble for species and the named ancestral nodes
values_traj <- values_traj[c(1:25,27,31,36,37,39),]

# Collate the clade info for each species so this can be appended to the tibble
clade_info <- c("Marsupialia", "Xenarthra", "Euarchontoglires", "Xenarthra", "Euarchontoglires",
                "Xenarthra", "Laurasiatheria", "Laurasiatheria", "Afrotheria", "Laurasiatheria",
                "Euarchontoglires", "Marsupialia", "Euarchontoglires", "Monotremata",
                "Laurasiatheria", "Marsupialia", "Euarchontoglires", "Afrotheria", "Marsupialia",
                "Marsupialia", "Laurasiatheria", "Marsupialia", "Ancestral", 
                "Ancestral", "Ancestral", "Ancestral", "Ancestral",
                "Ancestral", "Ancestral", "Ancestral")

# Collate all the species info so this can be appended to the tibble
species_info <- c(species_list, "Ancestral mammal", "Ancestral therian mammal", "Ancestral placental mammal", 
                "Ancestral Laurasiatheria", "Ancestral Euarchontoglires", 
                "Ancestral Afrotheria", "Ancestral Xenarthra", "Ancestral Marsupialia")

# Collate all the extant/ancestral info so this can be appended to the tibble
extant_info <- c("Extant","Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Extant",
                 "Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Extant",
                 "Extant", "Extant", "Extant", "Extant", "Extant", "Extant", "Ancestral", "Ancestral", 
                 "Ancestral", "Ancestral", "Ancestral", "Ancestral", "Ancestral", "Ancestral")

# Append the specimen info to the tibble
values_traj <- values_traj %>% as_tibble() %>% mutate(Species_Phylogeny = rownames(values_traj), Species = species_info, Clade = clade_info, Extant = extant_info)

# Order the data by clade and convert to factor for plotting and plot legend
values_traj$Clade <- as.factor(values_traj$Clade)
values_traj <- values_traj %>% arrange(Clade)


# Create a species colour palette
mypalette_species <- c("mediumpurple1", "darkorchid4", # Afrotheria
                       "gray10", "gray35", "gray55", "darkorange2", "green3", "mediumpurple", "palevioletred", "dodgerblue2", # Ancestral
                       "darkolivegreen1", "lawngreen", "green3", "chartreuse4", "darkgreen", # Euarchontoglires
                       "burlywood1", "sandybrown", "darkorange2", "orangered1", "red3", # Laurasiatheria
                       "paleturquoise1", "turquoise1", "cyan3", "dodgerblue2", "blue2", "midnightblue", # Marsupialia
                       "gold1", # Monotremata
                       "pink1", "palevioletred", "mediumvioletred") # Xenarthra
image(1:28, 1, as.matrix(1:28), col = mypalette_species, xlab = "Macroscelides proboscideus, Setifer setosus, Ancestral mammal, Ancestral Laurasiatheria, Ancestral Euarchontoglires, Ancestral Afrotheria, Ancestral Xenarthra, Ancestral Marsupialia, Cebus apella, Dasyprocta leporina, Microcebus murinus, Mus musculus, Rattus rattus, Epomops franqueti, Felis catus, Manis tricuspis, Phacochoerus aethiopicus, Talpa europaea, Bettongia penicillata, Monodelphis domestica, Phascolarctos cincereus, Setonix brachyurus, Sminthopsis macroura, Trichosaurus vulpecula, Ornithorhynchus anatinus, Bradypus tridactylus, Cyclopes didactylus, Dasypus novemcinctus",
      ylab = "", yaxt = "n")

# Create a clade colour palette
mypalette_species_clade <- c("mediumpurple", "mediumpurple", # Afrotheria
                       "gray10", "gray35", "gray55", "darkorange2", "green3", "mediumpurple", "palevioletred", "dodgerblue2", # Ancestral
                       "green3", "green3", "green3", "green3", "green3", # Euarchontoglires
                       "darkorange2", "darkorange2", "darkorange2", "darkorange2", "darkorange2", # Laurasiatheria
                       "dodgerblue2", "dodgerblue2", "dodgerblue2", "dodgerblue2", "dodgerblue2", "dodgerblue2", # Marsupialia
                       "gold1", # Monotremata
                       "palevioletred", "palevioletred", "palevioletred") # Xenarthra
image(1:28, 1, as.matrix(1:28), col = mypalette_species_clade, xlab = "Macroscelides proboscideus, Setifer setosus, Ancestral mammal, Ancestral Laurasiatheria, Ancestral Euarchontoglires, Ancestral Afrotheria, Ancestral Xenarthra, Ancestral Marsupialia, Cebus apella, Dasyprocta leporina, Microcebus murinus, Mus musculus, Rattus rattus, Epomops franqueti, Felis catus, Manis tricuspis, Phacochoerus aethiopicus, Talpa europaea, Bettongia penicillata, Monodelphis domestica, Phascolarctos cincereus, Setonix brachyurus, Sminthopsis macroura, Trichosaurus vulpecula, Ornithorhynchus anatinus, Bradypus tridactylus, Cyclopes didactylus, Dasypus novemcinctus",
      ylab = "", yaxt = "n")

# Order the dataset to match the species colours
values_traj$Species <- factor(values_traj$Species,
                             levels = c("Macroscelides proboscideus", "Setifer setosus", "Ancestral mammal",
                              "Ancestral therian mammal", "Ancestral placental mammal",
                              "Ancestral Laurasiatheria", "Ancestral Euarchontoglires", "Ancestral Afrotheria",      
                              "Ancestral Xenarthra", "Ancestral Marsupialia", "Cebus apella",            
                              "Dasyprocta leporina", "Microcebus murinus", "Mus musculus",            
                              "Rattus rattus", "Epomops franqueti", "Felis catus",       
                              "Manis tricuspis", "Phacochoerus aethiopicus", "Talpa europaea",          
                              "Bettongia penicillata", "Monodelphis domestica", "Phascolarctos cincereus", 
                              "Setonix brachyurus", "Sminthopsis macroura", "Trichosaurus vulpecula",   
                              "Ornithorhynchus anatinus", "Bradypus tridactylus", "Cyclopes didactylus",      
                              "Dasypus novemcinctus"),
                               ordered = TRUE)


############################################################################################################

# STEP 7: Plot the species and ancestral allometric trajectories

# Plot allometric trajectories with abline for each group and ancestral node
allometry_anc_nodes_ggplot <- ggplot(allometry_plot_data, aes(x = logCsize, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = values_traj, aes(intercept = Intercept, slope = Slope, colour = Species, linetype = Extant, alpha = Extant), size = 1)+
  scale_colour_manual(values = mypalette_species_clade)+    
  scale_alpha_discrete(range = c(1, 0.2))+
  scale_linetype_manual(name = NULL, labels = c("Ancestral", "Extant"), values = c(2,1))+
  theme_classic(base_size = 12)+
  ggtitle("Species allometry and ancestral states allometery")+
  xlab("logCS")+
  ylab("RegScores")+
  ylim(-1.2,1.2)+
  xlim(2,8)
allometry_anc_nodes_ggplot


#####

# Plot the anc.ML traits on the phylogeny

# Import trees in Nexus format - branch lengths needed
tree <- "Data/my_mammal_tree.nexus"  
# Read the tree for analysis
tree_species <- read.nexus(tree) 
Phylogeny_species_list <- read.csv("Data/tree_taxa_names.csv")

values_traj_slope <- as.matrix(values_traj$Slope[1:22])
names(values_traj_slope) <- Phylogeny_species_list$Taxa_names

values_traj_intercept <- as.matrix(values_traj$Intercept[1:22])
names(values_traj_intercept) <- Phylogeny_species_list$Taxa_names

# Phylogeny map for the slope
obj <- contMap(tree_species, values_traj_slope, plot=FALSE)
plot(obj, legend = 0.7 * max(nodeHeights(tree_species)), fsize=c(0.7,0.9))

# Phylogeny map for the intercept
obj <- contMap(tree_species, values_traj_intercept, plot=FALSE)
plot(obj, legend = 0.7 * max(nodeHeights(tree_species)), fsize=c(0.7,0.9))


### END ###
