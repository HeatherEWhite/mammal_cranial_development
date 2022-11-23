# Title: Variably present bones

# Name: Heather White

# Date created: 20/08/21

# Last modified: 18/11/22

# License: MIT license


# Dealing with LMs for bones that are variably present
# Fixing LM position into one place for all landmarks associated with the same variably present bone

#######################################################################################################


rm(list = ls())

library(geomorph)
library(Morpho)
library(rgl)


#######################################################################################################

# STEP 1: Load the data

# Load the mirrored, Procrusted, and RHS deleted LM data (LHS)
load(file = "Data/LHS_only_Procrusted_LMs_mirrored.Rdata")
# Missing LMs (9999) already dealt with here 
# Only need to deal with variably present bones


####################################################################################################

# STEP 2: Dealing with ABSENT bones (variably present) after Procrustes 

# This approach slide LMs back on top of each other

# Create a new array containing the Procrusted LMs, to later manipulate the absent bone LMs
absent_Proc <- final_LM_data

# Read in a .csv containing the absent bone information
# In this .csv you only need to include the first LM for each absent bone
absent<-read.csv("Data/absent_LMs.csv")
# Read in a .csv containing details on which LMs are associated with which bones
LM_table <- read.csv("Data/LM_list.csv")

# Assign the LMs for the missing bones to variables
lm_jugal <- LM_table$LM[which(LM_table$Bone%in%c("jugal"))]
lm_ventral_premax <- LM_table$LM[which(LM_table$Bone%in%c("ventral_premax"))]
lm_IP <- LM_table$LM[which(LM_table$Bone%in%c("interparietal"))]
lm_paraoccipital <- LM_table$LM[which(LM_table$Bone%in%c("paraoccipital"))]
lm_premax <- LM_table$LM[which(LM_table$Bone%in%c("premax"))]
lm_basisphenoid <- LM_table$LM[which(LM_table$Bone%in%c("basisphenoid"))]


# Amend landmarks for the variably present bones

# In the variably present jugal - make LMs 16-19 all the same for specimens missing the jugal
for (i in 1:nrow(absent)){
  if( !is.na(absent$Jugal[i]))
    absent_Proc[lm_jugal,c(1:3),i] <- matrix(final_LM_data[16,c(1:3),i], nrow = length(lm_jugal), ncol=3,byrow=TRUE)
}
# In the variably present ventral premax - make LMs 67-69 all the same for specimens missing the ventral premax
for (i in 1:nrow(absent)){
  if( !is.na(absent$Ventral_premax[i]))
    absent_Proc[lm_ventral_premax,c(1:3),i] <- matrix(final_LM_data[67,c(1:3),i], nrow = length(lm_ventral_premax), ncol=3,byrow=TRUE)
}
# In the variably present interparietal- make LMs 38-40 all the same for specimens missing the interparietal
for (i in 1:nrow(absent)){
  if( !is.na(absent$Interparietal[i]))
    absent_Proc[lm_IP,c(1:3),i] <- matrix(final_LM_data[38,c(1:3),i], nrow = length(lm_IP), ncol=3,byrow=TRUE)
}
# In the variably present paraoccipital - make LMs 46-47 all the same for specimens missing the paraoccipital
for (i in 1:nrow(absent)){
  if( !is.na(absent$Paraoccipital[i]))
    absent_Proc[lm_paraoccipital,c(1:3),i] <- matrix(final_LM_data[46,c(1:3),i], nrow = length(lm_paraoccipital), ncol=3,byrow=TRUE)
}
# In the variably present premax - make LMs 5-8 all the same for specimens missing the premax
for (i in 1:nrow(absent)){
  if( !is.na(absent$Premaxilla[i]))
    absent_Proc[lm_premax,c(1:3),i] <- matrix(final_LM_data[5,c(1:3),i], nrow = length(lm_premax), ncol=3,byrow=TRUE)
}
# In the variably present basisphenoid - make LMs 55-57 all the same for specimens missing the basisphenoid
for (i in 1:nrow(absent)){
  if( !is.na(absent$Basisphenoid[i]))
    absent_Proc[lm_basisphenoid,c(1:3),i] <- matrix(final_LM_data[55,c(1:3),i], nrow = length(lm_basisphenoid), ncol=3,byrow=TRUE)
}


# Check the adjusted landmarks
# These should now all be identical for the landmarks associated with variably present bones

# Check all the missing LMs are now the same for the jugal
absent_Proc[,,24]
absent_Proc[lm_jugal,,24]
# Check for the ventral_premax
absent_Proc[,,42]
absent_Proc[lm_ventral_premax,,42]
# Check for the interparietal
absent_Proc[,,11]
absent_Proc[lm_IP,,11]
# Check for the paraoccipital
absent_Proc[,,13]
absent_Proc[lm_paraoccipital,,13]
# Check for the premax
absent_Proc[,,40]
absent_Proc[lm_premax,,40]
# Check for the basisphenoid
absent_Proc[,,77]
absent_Proc[lm_basisphenoid,,77]


save(absent_Proc, file = "Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")


### END ###
