# Title: Mirroring LMs

# Name: Heather White

# Date created:09/08/21

# Last modified: 18/11/22

# License: MIT license


# Mirroring cranial left-side landmarks across the midline to the right-side
# Mirroring performed using 3 midline LMs only
# Code adapted from Andrew Knapp and Sandra Álvarez-Carretero

#######################################################################################################

rm(list = ls())

library(geomorph)
library(Morpho)
library(rgl)


#######################################################################################################

# STEP 1: Load the data

# Locate the folder containing the cranial landmarks
folder <- 'Data/LMs'


# Read in the cranial landmark data into an array
nspecimens<-165 # number of specimens within the dataset
# Create a list containing specimen names - from the .pts files
ptslist<-dir(path = folder, pattern='*.pts',recursive=T)

# Set the working directory to the LM folder - not the overall R project
setwd("/Users/heatherwhite/Documents/PhD/PhD/R_Projects/mammal_cranial_development/Data/LMs")

# Create and fill an array containing the .pts landmarks
ptsarray<-array(dim=c(69,3,165)) # dim=c(number of landmarks on the LHS (to mirror), number of dimensions, number of specimens)
for(i in 1:length(ptslist))
{
  ptsarray[,,i]<-as.matrix(read.table(file=ptslist[i],skip=2,header=F,sep="",row.names=1))
}


# Return to the R project directory
setwd("/Users/heatherwhite/Documents/PhD/PhD/R_Projects/mammal_cranial_development")

# Label the array with the specimen name for each matrix within the array 
species_names <- as.data.frame(ptslist)
dimnames(ptsarray)[3] = species_names
lm_array<-ptsarray



#######################################################################################################

# STEP 2: Deal with the missing LMs


# Estimate the value of the missing LMs
lm_array[which(lm_array==9999)] <- NA
lm_array <- estimate.missing(lm_array,method="TPS")


#######################################################################################################

# STEP 3: Define matrix for vectors to get plane to reflect

# This creates a matrix with 3 positions for each specimen 
# This matrix will then be filled by specimen specific midline landmarks

# Create an empty 165 x 3 matrix (number of specimens x midline landmark numbers)
mat.vectors.plane <- matrix( rep(0, length(dimnames(lm_array)[[3]])*3), nrow = length(dimnames(lm_array)[[3]]), ncol = 3 )
rownames( mat.vectors.plane ) <- as.character( dimnames( lm_array )[[3]] )
colnames( mat.vectors.plane ) <- c( "v1", "v2", "v3" )


#######################################################################################################

# STEP 4: Append midline LMs to the mat.vectors.plane matrix


# These three numbers are the 3 landmarks that will designate the skull midline #####
v.all.specimens <- c(41,51,54)

# Apply midline landmarks to all specimens
for ( i in seq(1:length(dimnames(lm_array)[[3]])) ){
  mat.vectors.plane[i,] <- v.all.specimens
}


# Add exceptions for groups of specimens (if neccessary)
Dasyprocta <- c(41, 51, 59)
Manis_Microcebus_Mus_Talpa <- c(41, 51, 67)
Monodelphis_Phascolarctos_Trichosaurus <- c(41, 51, 52)
Ornithorhynchus <- c(25, 51, 65)
for ( i in c(28:35) ){
mat.vectors.plane[i,] <- Dasyprocta
}
for ( i in c(60:66, 67:70, 84:89, 149:156) ){
  mat.vectors.plane[i,] <- Manis_Microcebus_Mus_Talpa
}
for ( i in c(71:83, 105:115, 157:165) ){
  mat.vectors.plane[i,] <- Monodelphis_Phascolarctos_Trichosaurus
}
for ( i in c(90:93) ){
  mat.vectors.plane[i,] <- Ornithorhynchus
}


# Add exceptions for individual specimens (if neccessary)
Bettongia_adult <- c(20, 51, 68)
mat.vectors.plane[10,] <- Bettongia_adult
Bradypus_1867.4.12.579 <- c(20, 43, 59)
mat.vectors.plane[11,] <- Bradypus_1867.4.12.579
Bradypus_1989.226 <- c(5, 41, 52)
mat.vectors.plane[12,] <- Bradypus_1989.226
Bradypus_51.1173 <- c(4, 42, 51)
mat.vectors.plane[14,] <- Bradypus_51.1173
Bradypus_no_number <- c(5, 42, 57)
mat.vectors.plane[15,] <- Bradypus_no_number
Bradypus_9195 <- c(42, 51, 68)
mat.vectors.plane[17,] <- Bradypus_9195
Cebus_12.6.5.8 <- c(4, 20, 42)
mat.vectors.plane[18,] <- Cebus_12.6.5.8
Cebus_28.2.9.3 <- c(1, 20, 51)
mat.vectors.plane[19,] <- Cebus_28.2.9.3
Cebus_67.4.12.394 <- c(42, 51, 54)
mat.vectors.plane[21,] <- Cebus_67.4.12.394
Cebus_71.3174 <- c(4, 42, 52)
mat.vectors.plane[22,] <- Cebus_71.3174
Cebus_71.3177 <- c(1, 20, 43)
mat.vectors.plane[23,] <- Cebus_71.3177
Cyclopes_2010.105 <- c(5, 20, 51)
mat.vectors.plane[24,] <- Cyclopes_2010.105
Cyclopes_55.1.26.345 <- c(42, 51, 65)
mat.vectors.plane[26,] <- Cyclopes_55.1.26.345
Cyclopes_26.12.4.68 <- c(19, 20, 51)
mat.vectors.plane[25,] <- Cyclopes_26.12.4.68
Cyclopes_6XII1906a <- c(25, 51, 64)
mat.vectors.plane[27,] <- Cyclopes_6XII1906a
Dasyprocta_12.5.11.8 <- c(5, 20, 57)
mat.vectors.plane[28,] <- Dasyprocta_12.5.11.8
Dasyprocta_67.4.12.472 <- c(38, 51, 59)
mat.vectors.plane[30,] <- Dasyprocta_67.4.12.472
Dasyprocta_67.4.12.497 <- c(38, 51, 68)
mat.vectors.plane[31,] <- Dasyprocta_67.4.12.497
Dasyprocta_67.4.12.585_spec2 <- c(38, 51, 67)
mat.vectors.plane[33,] <- Dasyprocta_67.4.12.585_spec2
Dasyprocta_97.8.13.2 <- c(20, 51, 67)
mat.vectors.plane[34,] <- Dasyprocta_97.8.13.2
Dasyprocta_2021.1 <- c(38, 51, 58)
mat.vectors.plane[35,] <- Dasyprocta_2021.1
Dasypus_1901.331 <- c(42, 51, 68)
mat.vectors.plane[36,] <- Dasypus_1901.331
Dasypus_1910.409b <- c(25, 51, 65)
mat.vectors.plane[37,] <- Dasypus_1910.409b
Dasypus_1910.410 <- c(25, 51, 68)
mat.vectors.plane[38,] <- Dasypus_1910.410
Dasypus_z.134 <- c(20, 58, 64)
mat.vectors.plane[39,] <- Dasypus_z.134
Dasypus_40651 <- c(4, 42, 51)
mat.vectors.plane[41,] <- Dasypus_40651
Epomops_66.3499 <- c(1, 26, 51)
mat.vectors.plane[43,] <- Epomops_66.3499
Epomops_68.353 <- c(4, 38, 51)
mat.vectors.plane[44,] <- Epomops_68.353
Epomops_80.7.21.4 <- c(1, 20, 57)
mat.vectors.plane[46,] <- Epomops_80.7.21.4
Felis_19.7.7.3514 <- c(20, 51, 68)
mat.vectors.plane[47,] <- Felis_19.7.7.3514
Felis_1952.10.20.1 <- c(20, 51, 68)
mat.vectors.plane[48,] <- Felis_1952.10.20.1
Felis_1952.10.20.2 <- c(1, 20, 51)
mat.vectors.plane[49,] <- Felis_1952.10.20.2
Felis_1992.184 <- c(38, 51, 54)
mat.vectors.plane[51,] <- Felis_1992.184
Macroscelides_2.9.1.17 <- c(1, 51, 54)
mat.vectors.plane[54,] <- Macroscelides_2.9.1.17
Macroscelides_2.9.1.18 <- c(20, 51, 59)
mat.vectors.plane[55,] <- Macroscelides_2.9.1.18
Macroscelides_48 <- c(38, 51, 59)
mat.vectors.plane[57,] <- Macroscelides_48
Macroscelides_w3 <- c(20, 51, 65)
mat.vectors.plane[58,] <- Macroscelides_w3
Macroscelides_w39 <- c(4, 38, 54)
mat.vectors.plane[59,] <- Macroscelides_w39
Manis_1.8.9.108 <- c(4, 51, 54)
mat.vectors.plane[60,] <- Manis_1.8.9.108
Manis_1999.102 <- c(38, 51, 65)
mat.vectors.plane[61,] <- Manis_1999.102
Manis_1999.93 <- c(38, 51, 54)
mat.vectors.plane[62,] <- Manis_1999.93
Manis_66.3562 <- c(38, 54, 67)
mat.vectors.plane[63,] <- Manis_66.3562
Manis_9.1.4.66 <- c(38, 51, 57)
mat.vectors.plane[64,] <- Manis_9.1.4.66
Microcebus_7030 <- c(20, 51, 65)
mat.vectors.plane[67,] <- Microcebus_7030
Microcebus_infant1 <- c(38, 52, 67)
mat.vectors.plane[68,] <- Microcebus_infant1
Microcebus_infant2 <- c(19, 51, 68)
mat.vectors.plane[69,] <- Microcebus_infant2
Monodelphis_20days <- c(4, 51, 68)
mat.vectors.plane[72,] <- Monodelphis_20days
Monodelphis_22days <- c(4, 51, 67)
mat.vectors.plane[73,] <- Monodelphis_22days
Monodelphis_25days <- c(4, 57, 67)
mat.vectors.plane[74,] <- Monodelphis_25days
Monodelphis_30days <- c(4, 51, 65)
mat.vectors.plane[75,] <- Monodelphis_30days
Monodelphis_35days <- c(4, 51, 59)
mat.vectors.plane[76,] <- Monodelphis_35days
Monodelphis_57days <- c(4, 51, 67)
mat.vectors.plane[78,] <- Monodelphis_57days
Mus_P14 <- c(38, 51, 67)
mat.vectors.plane[84,] <- Mus_P14
Mus_P18 <- c(1, 38, 51)
mat.vectors.plane[85,] <- Mus_P18
Mus_P25 <- c(38, 51, 68)
mat.vectors.plane[86,] <- Mus_P25
Mus_4.5months <- c(19, 51, 67)
mat.vectors.plane[87,] <- Mus_4.5months
Mus_5wk <- c(19, 51, 67)
mat.vectors.plane[88,] <- Mus_5wk
Mus_9wk <- c(38, 51, 68)
mat.vectors.plane[89,] <- Mus_9wk
Ornithorhynchus_no_number <- c(26, 51, 65)
mat.vectors.plane[93,] <- Ornithorhynchus_no_number
Phacochoerus_66.428<- c(20, 51, 64)
mat.vectors.plane[95,] <- Phacochoerus_66.428
Phacochoerus_66.523 <- c(20, 51, 64)
mat.vectors.plane[97,] <- Phacochoerus_66.523
Phacochoerus_71.2125 <- c(20, 51, 59)
mat.vectors.plane[103,] <- Phacochoerus_71.2125
Phacochoerus_71.7.3.4 <- c(1, 41, 51)
mat.vectors.plane[104,] <- Phacochoerus_71.7.3.4
Phascolarctos_57days <- c(38, 51, 54)
mat.vectors.plane[106,] <- Phascolarctos_57days
Phascolarctos_51days <- c(38, 51, 54)
mat.vectors.plane[107,] <- Phascolarctos_51days
Phascolarctos_PcinLge02 <- c(4, 42, 54)
mat.vectors.plane[113,] <- Phascolarctos_PcinLge02
Phascolarctos_PcinLge03 <- c(20, 51, 67)
mat.vectors.plane[114,] <- Phascolarctos_PcinLge03
Phascolarctos_PcinLge05 <- c(38, 51, 58)
mat.vectors.plane[115,] <- Phascolarctos_PcinLge05
Rattus_1997.69 <- c(1, 38, 51)
mat.vectors.plane[116,] <- Rattus_1997.69
Rattus_1999.14 <- c(4, 42, 54)
mat.vectors.plane[117,] <- Rattus_1999.14
Rattus_52.1111 <- c(19, 42, 54)
mat.vectors.plane[118,] <- Rattus_52.1111
Rattus_70.103 <- c(51, 64, 68)
mat.vectors.plane[119,] <- Rattus_70.103
Rattus_79.1343 <- c(19, 43, 65)
mat.vectors.plane[120,] <- Rattus_79.1343
Setifer_1974.484 <- c(1, 41, 54)
mat.vectors.plane[121,] <- Setifer_1974.484
Setifer_55.12.26.304 <- c(4, 51, 68)
mat.vectors.plane[122,] <- Setifer_55.12.26.304
Setifer_70.360 <- c(20, 51, 67)
mat.vectors.plane[123,] <- Setifer_70.360
Setifer_74.545 <- c(4, 51, 67)
mat.vectors.plane[124,] <- Setifer_74.545
Setifer_74.554 <- c(42, 51, 65)
mat.vectors.plane[125,] <- Setifer_74.554
Setifer_76.273 <- c(20, 51, 67)
mat.vectors.plane[126,] <- Setifer_76.273
Setifer_76.274 <- c(1, 42, 54)
mat.vectors.plane[127,] <- Setifer_76.274
Setifer_79.545 <- c(4, 41, 54)
mat.vectors.plane[128,] <- Setifer_79.545
Setonix_1989.585 <- c(38, 52, 67)
mat.vectors.plane[134,] <- Setonix_1989.585
Setonix_1989.587 <- c(19, 51, 68)
mat.vectors.plane[136,] <- Setonix_1989.587
Setonix_1989.590 <- c(4, 38, 67)
mat.vectors.plane[139,] <- Setonix_1989.590
Setonix_6.8.1.245 <- c(19, 41, 51)
mat.vectors.plane[141,] <- Setonix_6.8.1.245
Sminthopsis_64days <- c(42, 51, 59)
mat.vectors.plane[143,] <- Sminthopsis_64days
Sminthopsis_adult <- c(19, 51, 67)
mat.vectors.plane[144,] <- Sminthopsis_adult
Sminthopsis_19days <- c(42, 51, 67)
mat.vectors.plane[145,] <- Sminthopsis_19days
Sminthopsis_74days <- c(4, 51, 68)
mat.vectors.plane[148,] <- Sminthopsis_74days
Talpa_19.4.57 <- c(42, 51, 68)
mat.vectors.plane[149,] <- Talpa_19.4.57
Talpa_57.367 <- c(4, 51, 67)
mat.vectors.plane[152,] <- Talpa_57.367
Talpa_57.368 <- c(51, 54, 68)
mat.vectors.plane[153,] <- Talpa_57.368
Talpa_57.370 <- c(1, 42, 51)
mat.vectors.plane[154,] <- Talpa_57.370
Talpa_57.371 <- c(4, 41, 68)
mat.vectors.plane[155,] <- Talpa_57.371
Trichosaurus_TV1 <- c(4, 51, 67)
mat.vectors.plane[157,] <- Trichosaurus_TV1
Trichosaurus_TV2 <- c(26, 51, 67)
mat.vectors.plane[158,] <- Trichosaurus_TV2
Trichosaurus_TV3 <- c(1, 38, 57)
mat.vectors.plane[159,] <- Trichosaurus_TV3
Trichosaurus_TV8 <- c(38, 51, 65)
mat.vectors.plane[162,] <- Trichosaurus_TV8
Trichosaurus_15days <- c(38, 51, 54)
mat.vectors.plane[164,] <- Trichosaurus_15days
Trichosaurus_21days <- c(38, 51, 54)
mat.vectors.plane[165,] <- Trichosaurus_21days
Manis_95.7.17.1 <- c(20, 51, 65)
mat.vectors.plane[66,] <- Manis_95.7.17.1
Phacochoerus_71.2124 <- c(20, 43, 58)
mat.vectors.plane[102,] <- Phacochoerus_71.2124
Ornithorhynchus_1859.6.30.13 <- c(20, 43, 65)
mat.vectors.plane[90,] <- Ornithorhynchus_1859.6.30.13


# Visualise the midline landmark matrix
mat.vectors.plane

# Check which specimen is associated which each row
row.names(mat.vectors.plane)


#######################################################################################################

# STEP 5 - Implement the function to reflect the landmarks
# All you need to do is run this function - don't change any inputs etc.

#---------------------------------------------------------------------------------------#
# FUNCTION reflect.lmks (Now it is also saved in Functions.R)                           #
#=======================================================================================#
#
# X = Array of size "k x q x n", where "k" is the number
#     landmarks, "q" the number of coordinates, and "n" the
#     number of specimens.
#     Note that this array should not have NA values.
#     If it has NA values, please run before:
#        estimate.missing( A , method = c("TPS", "Reg" ) )
# v = Matrix of length "n x 3", where "n" is the number of 
#     specimens and "3" the number of points to create a plane 
#     thanks to which the landmarks are reflected in 
#     Morpho::mirror2plane. If the three points are the same 
#     for all specimens (each row same values), you can just add
#     a vector of length 3 with these points.
# midline = Vector with the position of the coordinates that
#           are part of the midline
# plot.res = Boolean, TRUE if you want plots, FALSE otherwise.
#            NOTE: You will have as many 3D plots as specimens,
#            so be aware of that before using T

reflect.lmks <- function( X, v, midline, plot.res ){ 
  
  #\\ Get number of specimens and position 
  #\\ of lmks
  sp <- dim( X )[3]
  if( length( v ) == 3 ){
    v1 <- v[1]
    v2 <- v[2]
    v3 <- v[3]
  } else {
    v1 <- v[,1]
    v2 <- v[,2]
    v3 <- v[,3]
  }
  
  #\\ Get reflected landmarks
  # 1. Create array in which final reflected
  #    landmarks will be saved. Note that midlines
  #    are not added as the mean will be used
  reflect.out <- array( dim = c( lmks   = dim( X )[1],
                                 coords = dim( X )[2],
                                 sp     = sp ) )
  dimnames( reflect.out ) <- list( paste( "lmk.ref.",
                                          seq( 1:(dim( X )[1]) ), sep="" ),
                                   c( "x", "y", "z" ),
                                   c( unlist( dimnames( X )[3] ) ) )
  # 2. Create array in which mean midline will
  #    replace original midline to aid visualization
  X.mean <- X
  # 3. Get the reflected landmarks for all "i" specimens
  #    and fill in "reflect.out"
  for( i in seq( 1:sp ) ){
    
    #\\ Run Morpho::mirror2plane function for specimen "i"
    if( length( v ) == 3 ){
      tmp.new.mat <- Morpho::mirror2plane( x = X[,,i],
                                           v1 = X[v1,,i], v2 = X[v2,,i], v3 = X[v3,,i] )
    } else{ 
      tmp.new.mat <- Morpho::mirror2plane( x = X[,,i],
                                           v1 = X[v1[i],,i], v2 = X[v2[i],,i],
                                           v3 = X[v3[i],,i] )
    }
    
    #\\ Get mean midline to add it to "X.mean" for sp "i"
    x.mean.new.mat <- apply( cbind( X[midline,1,i], tmp.new.mat[midline,1]), 1, mean )
    y.mean.new.mat <- apply( cbind( X[midline,2,i], tmp.new.mat[midline,2]), 1, mean )
    z.mean.new.mat <- apply( cbind( X[midline,3,i], tmp.new.mat[midline,3]), 1, mean )
    mean.midline <- cbind( x.mean.new.mat, y.mean.new.mat, z.mean.new.mat )
    
    #\\ Change midline in "X.mean" by the mean of the midline
    #\\ previously calculated. Also remove midline coords in tmp.new.mat
    #\\ to add to final "reflect.out" array
    X.mean[midline,,i] <- mean.midline
    #reflect.out[,,i]   <- tmp.new.mat[-midline,] # done before, now I don't
    # to keep lmk.ref.X labels
    reflect.out[,,i]   <- tmp.new.mat
  }
  
  #\\ Remove midline coordinates from reflect out
  reflect.out <- reflect.out[-midline,,]
  
  #\\ Combine both data sets so we can plot them in 3D
  # 1. Create array for combined data set with original 
  #    midline and fill in
  comb.lmks <- dim( X )[1]*2 - length( midline )
  comb.arr  <- array( dim = c( lmks = comb.lmks, coords = 3, sp = sp ) )
  dimnames( comb.arr ) <- list( c( paste( "lmk",
                                          seq( 1:(dim( X )[1]) ), sep="" ),
                                   unlist( dimnames( reflect.out )[1] ) ),
                                c( "x", "y", "z" ),
                                c( unlist( dimnames( X )[3] ) ) )
  # 2. Create array for combined data set with mean 
  #    midline and fill in
  comb.arr.mean  <- comb.arr
  # 3. Fill in arrays
  comb.arr[1:(dim( X )[1]),,]             <- X
  comb.arr[((dim( X )[1])+1):comb.lmks,,] <- reflect.out
  comb.arr.mean[1:(dim( X )[1]),,]             <- X.mean
  comb.arr.mean[((dim( X )[1])+1):comb.lmks,,] <- reflect.out
  
  #\\ Plot 3D if user sets plot.res = T
  if( plot.res == T ){
    # Define colours so the mean midline of is coloured in
    # blue, the digitized side in black and the reflect side 
    # in dark green
    colours <- rep( "black", comb.lmks )
    colours[ midline ] <- "blue"
    colours[ v ] <- "red"
    colours[ dim( X )[1]:comb.lmks ] <- "darkgreen"
    
    # Plot combined data set
    for( i in seq( 1:sp ) ){
      open3d()
      plot3d( comb.arr.mean[,,i], col = colours, size = 10, aspect = FALSE,
              main = unlist( dimnames( X )[3] )[i] )
    }
  }
  
  #\\ Return both arrays
  return( list( original = comb.arr, mean = comb.arr.mean ) )
  
}
#=======================================================================================#
# END OF FUNCTION reflect.lmks                                                          #
#---------------------------------------------------------------------------------------#

#######################################################################################################

# STEP 6 - Get the arguments for function "reflect.lmks" and run function

# Assign the imported landmarks to variable X
X <- lm_array
# Assign the midline landmarks to use for the dataset to variable v
v <- mat.vectors.plane 

# List ALL landmarks that lie on the midline in EVERY specimen
# Therefore the landmarks that are in mat.vectors.plane might not all feature in here - if they aren't on the midline in every specimen
midline <- c(41,42,43,51,52)

# Implement the reflect.lmks function
comb.dataset <- reflect.lmks(X = X, v = v, midline = midline, plot.res = F) # Mark as T/F if do/don't want all the visual windows to open


#### Check the results

# Visualise the results on a specimen
# The below line has to be for the ASCII .ply file
scan2=ply2mesh(file="Data/plys_ascii/HW_Manis_tricuspis_91.363_skull_decimated.ply") 
# Check which number this spec is using the ptslist variable - have to input this number below, here as 65

# Highlight all four lines and run togethers to plot LMs on the skull
shade3d(scan2,col='white', specula = "black")
# To visualise both LHS, RHS and midline LMs on specimen
spheres3d(comb.dataset$original[c(1:40,44:50,53:69),,65],col=5,radius=0.2) #LHS landmarks
spheres3d(comb.dataset$original[c(70:133),,65],col=6,radius=0.2) # RHS landmarks
spheres3d(comb.dataset$original[c(41,42,43,51,52),,65],col=4,radius=0.2) #Midline landmarks



######################################################################################################

#  STEP 7: Perform Procrustes superimposition on the mirrored data and then delete the mirrored half

# Perform Procrustes
Procrusted_mirrored_skull <- gpagen(comb.dataset$original, verbose =T) 
# Subset out the coords
coords_data_mirrored <- Procrusted_mirrored_skull$coords 
# Subset out the centroid size
CS <- Procrusted_mirrored_skull$Csize 
# Access the mean shape data
meanshape <- mshape(Procrusted_mirrored_skull$coords) 
# Subset out the procrustes distances (Euclidean distance); in gpagen function verbose needs to be set to TRUE for this
ProcD <- Procrusted_mirrored_skull$procD 
# Convert Procrustes distances to a matrix
procdist <- as.matrix(ProcD) 

# Save the centroid size data
write.csv(CS, file = "Data/Centroid_size.csv")


# Checking if landmark ordering is correct
# Can't check these on a specimen, as they are the Procrusted LMs
# Can only visualise the LMs on the skull before Procrustes
spheres3d(coords_data_mirrored[1:69,,65], radius = 0.002, color = 'red')
spheres3d(coords_data_mirrored[70:133,,65], radius = 0.002, color = 'green')
spheres3d(coords_data_mirrored[midline,,65], radius = 0.002, color = 'black')


# After performing Procrustes remove the mirrored RHS
final_LM_data <- coords_data_mirrored[1:length(1:nrow(lm_array)),,]


# Save this final mirrored, Procrusted, and RHS deleted LM data
# To input this data into 'variably_present.R' code
save(final_LM_data, file = "Analysis/LHS_only_Procrusted_LMs_mirrored.Rdata")


### END ###
