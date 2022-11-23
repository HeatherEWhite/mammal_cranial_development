# Title: Plotting extreme PC morphologies

# Name: Heather White

# Date created: 25/10/21

# Last modified: 18/11/22

# License: MIT license


# Obtaining the LMs for the extreme PC shapes

#######################################################################################################

rm(list = ls())

library(Morpho)
library(geomorph)
library(rgl)

#######################################################################################################

# STEP 1: Load the data and format

# LMs mirrored (not Procrusted) - comb.dataset
load(file = "Data/Mirrored_LMs.Rdata")
# Getting the adults only - original data mirrored (no Procrustes)
adults_OG <- comb.dataset$original[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]
# Perform PCA for the non-Procrusted dataset
PCA1 <- gm.prcomp(adults_OG)


# LMs mirrored and Procrusted - coords_data_mirrored
load(file = "Data/mirrored_Procrusted_LMs.Rdata")
# Getting the adults only - mirrored and Procrusted data
adults <- coords_data_mirrored[,,c(10,11,18,25,29,39,43,48,56,60,67,82,87,91,94,113,117,121,141,144,150,157)]
# Perform PCA for the Procrustes mirrored dataset
PCA2 <- gm.prcomp(adults)
# PCA of full development dataset
PCA_all_both_sides <- gm.prcomp(coords_data_mirrored)


# Load the Procrusted coordinate data (LHS only) - only specimens
load(file ="Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")
# Perform PCA on Procrustes data
PCA_all <- gm.prcomp(absent_Proc)


# Load the Procrusted coordinate data (LHS only) - adults only
load(file = "Data/adult_Procrusted_LMs.Rdata")
# Perform PCA on Procrustes data
PCA_adults <- gm.prcomp(adults)


# Mesh to warp based on PC extremes
mesh.setonix <- ply2mesh(file = "Data/plys_ascii/Setonix_brachyurus_NHMUK_6.8.1.245.ply")


#######################################################################################################

# STEP 2: Obtain the warped meshes

# refmat = landmark configuration on a surface
# This is the mesh with the original mirrored data (not Procrusted otherwise the LMs don't sit on the mesh)

# tarmat = landmark configuration on a target surface
# Target configuation has to be the Procrusted data otherwise the extreme PCA shape will be incorrect

# Warps of the skull - adult dataset
max.mesh1 <- tps3d(mesh.setonix,adults_OG[,,19],PCA2$shapes$shapes.comp1$max,threads=1) # For PC1 max
shade3d(max.mesh1,col="white", specula ="black")
max.mesh1 <- tps3d(mesh.setonix,adults_OG[,,19],PCA2$shapes$shapes.comp1$min,threads=1) # For PC1 min
shade3d(max.mesh1,col="white", specula ="black")


# Warps of the skull - full dataset
max.mesh1 <- tps3d(mesh.setonix,comb.dataset$original[,,141],PCA_all_both_sides$shapes$shapes.comp1$max,threads=1) # For PC1 max
shade3d(max.mesh1,col="white", specula ="black")
max.mesh1 <- tps3d(mesh.setonix,comb.dataset$original[,,141],PCA_all_both_sides$shapes$shapes.comp1$min,threads=1) # For PC1 min
shade3d(max.mesh1,col="white", specula ="black")


#######################################################################################################

# STEP 3: Plot the extreme shapes - LMs only

# Plot the extreme shapes - full dataset
spheres3d(PCA_all$shapes$shapes.comp1$min, radius = 0.006, color = "grey") # Min shape for PC1
spheres3d(PCA_all$shapes$shapes.comp1$max, radius = 0.006, color = "grey") # Max shape for PC1 
spheres3d(PCA_all$shapes$shapes.comp2$min, radius = 0.006, color = "grey") # Min shape for PC2 
spheres3d(PCA_all$shapes$shapes.comp2$max, radius = 0.006, color = "grey") # Max shape for PC2

# Plot the extreme shapes - adults only
spheres3d(PCA_adults$shapes$shapes.comp1$min, radius = 0.006, color = "grey") # Min shape for PC1
spheres3d(PCA_adults$shapes$shapes.comp1$max, radius = 0.006, color = "grey") # Max shape for PC1 
spheres3d(PCA_adults$shapes$shapes.comp2$min, radius = 0.006, color = "grey") # Min shape for PC2 
spheres3d(PCA_adults$shapes$shapes.comp2$max, radius = 0.006, color = "grey") # Max shape for PC2


### END ###
