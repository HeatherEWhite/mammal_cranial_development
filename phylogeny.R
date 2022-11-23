# Creating a phylogeny

# Name: Heather White

# Date created: 01/10/21

# Last modified: 18/11/22

# License: MIT license


# Creating a phylogeny from Upham et al. 2019 for the 22 species included in the dataset of White et al. 2022

################################################################################################################

rm(list = ls())

library(ape)
library(geiger)


################################################################################################################

# STEP 1: Load the data

# Read in the full tree in nexus format from the literature
big_tree <- read.nexus("Data/Upham_NDexp_MCC.tre.txt")


################################################################################################################

# STEP 2: Produce the tree

# Keep the species within the dataset
# Names must match those within the tree (Upham et al. 2019)
my_tree<-keep.tip(big_tree, c("Ornithorhynchus_anatinus_ORNITHORHYNCHIDAE_MONOTREMATA", "Monodelphis_domestica_DIDELPHIDAE_DIDELPHIMORPHIA",
                                   "Sminthopsis_macroura_DASYURIDAE_DASYUROMORPHIA", "Phascolarctos_cinereus_PHASCOLARCTIDAE_DIPROTODONTIA",
                                   "Trichosurus_vulpecula_PHALANGERIDAE_DIPROTODONTIA", "Setonix_brachyurus_MACROPODIDAE_DIPROTODONTIA",
                                   "Bettongia_penicillata_POTOROIDAE_DIPROTODONTIA", "Dasypus_novemcinctus_DASYPODIDAE_CINGULATA", 
                                   "Bradypus_tridactylus_BRADYPODIDAE_PILOSA", "Cyclopes_didactylus_CYCLOPEDIDAE_PILOSA", 
                                   "Macroscelides_proboscideus_MACROSCELIDIDAE_MACROSCELIDEA", "Setifer_setosus_TENRECIDAE_AFROSORICIDA",
                                   "Talpa_europaea_TALPIDAE_EULIPOTYPHLA", "Epomops_franqueti_PTEROPODIDAE_CHIROPTERA", 
                                   "Phacochoerus_aethiopicus_SUIDAE_CETARTIODACTYLA", "Felis_catus_FELIDAE_CARNIVORA",
                                   "Phataginus_tricuspis_MANIDAE_PHOLIDOTA", "Rattus_rattus_MURIDAE_RODENTIA",
                                   "Microcebus_murinus_CHEIROGALEIDAE_PRIMATES", "Mus_musculus_MURIDAE_RODENTIA",
                                   "Dasyprocta_leporina_DASYPROCTIDAE_RODENTIA", "Sapajus_apella_CEBIDAE_PRIMATES"))

write.nexus(my_tree, file = 'Data/my_mammal_tree.nexus')


# Plot the timescaled phylogeny for this dataset
plot(my_tree, cex = 0.5)
axisPhylo()


######

# Produce a placental only tree
my_tree_placentals<-keep.tip(big_tree, c("Dasypus_novemcinctus_DASYPODIDAE_CINGULATA", 
                              "Bradypus_tridactylus_BRADYPODIDAE_PILOSA", "Cyclopes_didactylus_CYCLOPEDIDAE_PILOSA", 
                              "Macroscelides_proboscideus_MACROSCELIDIDAE_MACROSCELIDEA", "Setifer_setosus_TENRECIDAE_AFROSORICIDA",
                              "Talpa_europaea_TALPIDAE_EULIPOTYPHLA", "Epomops_franqueti_PTEROPODIDAE_CHIROPTERA", 
                              "Phacochoerus_aethiopicus_SUIDAE_CETARTIODACTYLA", "Felis_catus_FELIDAE_CARNIVORA",
                              "Phataginus_tricuspis_MANIDAE_PHOLIDOTA", "Rattus_rattus_MURIDAE_RODENTIA",
                              "Microcebus_murinus_CHEIROGALEIDAE_PRIMATES", "Mus_musculus_MURIDAE_RODENTIA",
                              "Dasyprocta_leporina_DASYPROCTIDAE_RODENTIA", "Sapajus_apella_CEBIDAE_PRIMATES"))

write.nexus(my_tree_placentals, file = 'Data/my_mammal_tree_placentals.nexus')


### END ###
