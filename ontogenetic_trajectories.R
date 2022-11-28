# Title: Ontogenetic allometric trajectories

# Name: Heather White

# Date created: 02/09/21

# Last modified: 18/11/22

# License: MIT license


# Calculating ontogenetic trajectories for each species using allometry regressions
# Adapted from Morris et al. (2019)

#######################################################################################################

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(geomorph)
library(abind)


#######################################################################################################

# STEP 1: Load the data 

# Load the PC scores for each specimen
load(file ="Data/PC_scores_all_specimens.Rdata")
rownames(PC_scores_PC1_flipped) <- PC_scores_PC1_flipped$Specimen
PC_scores_PC1_flipped <- PC_scores_PC1_flipped[,-165]

# Load the landmark coordinates and specimen info
load(file = "Data/LHS_absent_Procrusted_skull_LMs_slid.Rdata")
info <- read.csv(file = "Data/Specimen_info.csv", header =T, sep = ",")
# Remove specimen names from the coordinate dataset - needs to be as numbers for later analysis
absent_Proc <- unname(absent_Proc)


#######################################################################################################

# STEP 2: Format the data 

# Set colours
# Run PlottingValues_Function.R first
# This function from Morris et al. 2019
mammal.PlottingValues <- PlottingValues(X=info,ColorGroup = "Species",
                                        ShapeGroup = "Major_clades")


# Get the raw Procrustes data into a geomorph dataframe
gdf <- geomorph.data.frame(absent_Proc, CSize = info$CS, logCsize = info$CS_logged, 
                           Species = info$Species, Specimens = info$Specimen_name)
# name first part of data frame (containing Procrustes data)
names(gdf) <-c("Proc_coords","CSize", "logCsize", "Species", "Specimens")


# Assign colours to the species
col.species = c("paleturquoise1", "pink1", "darkolivegreen1", "palevioletred",
                "lawngreen", "mediumvioletred", "burlywood1", "sandybrown", "mediumpurple1",
                "darkorange2", "green3", "turquoise1", "chartreuse4", "gold1",
                "orangered1", "cyan3", "darkgreen", "darkorchid4",
                "dodgerblue2", "blue2", "red3", "midnightblue")[gdf$Species]


# Convert the mammal.PlottingValues to the colours I want - matching above
col.gp.1 <- c("paleturquoise1", "pink1", "darkolivegreen1", "palevioletred",
              "lawngreen", "mediumvioletred", "burlywood1", "sandybrown", "mediumpurple1",
              "darkorange2", "green3", "turquoise1", "chartreuse4", "gold1",
              "orangered1", "cyan3", "darkgreen", "darkorchid4",
              "dodgerblue2", "blue2", "red3", "midnightblue")
names( col.gp.1 ) <- levels(info$Species)
mammal.PlottingValues$legendcolor <- col.gp.1


# Perform PCA to get PCA results and output
PCA <- gm.prcomp(absent_Proc)
PCA_summary <- summary(PCA)
PCA_summary <- PCA_summary$PC.summary


# Subset dataset by species
GPAList=list()
CovariatesList=list()
SpeciesList=list()
CladeList=list()
PCList=list()
Shapes=list()
CSList=list()
LogCsize=list()
Colors=list()
Taxa=list()
for (g in 1:length(levels(info[["Species"]]))){
  x<-NULL
  x<-as.array(gdf$Proc_coords[,,grep(levels(info[["Species"]])[g],info[["Species"]])])
  y<-NULL
  y<-info[grep(levels(info[["Species"]])[g],info[["Species"]]),]
  z<-NULL
  z<-as.matrix(PC_scores_PC1_flipped[grep(levels(info[["Species"]])[g],info[["Species"]]),])
  {
    GPAList[[print(levels(info[["Species"]])[g])]]<-as.array(x)
    CSList[[levels(info[["Species"]])[g]]]<-info$CS[grep(levels(info[["Species"]])[g],info[["Species"]])]
    LogCsize[[levels(info[["Species"]])[g]]] <- info$CS_logged[grep(levels(info[["Species"]])[g],info[["Species"]])]
    CovariatesList[[levels(info[["Species"]])[g]]]<-y
    SpeciesList[[levels(info[["Species"]])[g]]]<-y$Species
    CladeList[[levels(info[["Species"]])[g]]]<-y$Major_clade
    Shapes[[levels(info[["Species"]])[g]]]<-mammal.PlottingValues$shape[grep(levels(info[["Species"]])[g],info[["Species"]])]
    Colors[[levels(info[["Species"]])[g]]]<-mammal.PlottingValues$color[grep(levels(info[["Species"]])[g],info[["Species"]])]
    Taxa<-names(GPAList)
  }
  if (ncol(z)==1) {
    z <- t(z)
    PCList[[levels(info[["Species"]])[g]]]<-z
  }else{
    PCList[[levels(info[["Species"]])[g]]]<-z
  }
}

# Remove rownames of PC matrix - otherwise this confuses later analysis
for (i in 1:22){
rownames(PCList[[i]]) <- NULL
}

# Subsetted data by species combined as one list
appended_info_species <- list(name = paste(info$Species[-23]), taxa = Taxa[-23], 
                              coords = GPAList[-23], PCvalues = PCList[-23], 
                              CSize = CSList[-23], logCsize = LogCsize[-23], 
                              covariates = CovariatesList[-23], species = SpeciesList[-23], 
                              clades = CladeList[-23], color = Colors[-23], shape = Shapes[-23])
A <- list()
for (i in 1:length(appended_info_species$taxa)){
  if (!is.na(dim(appended_info_species$coords[[i]])[3]) & dim(appended_info_species$coords[[i]])[3]>2) {
    A$PCvalues[[appended_info_species$taxa[i]]] <- appended_info_species$PCvalues[[i]]
    A$CSize[[appended_info_species$taxa[i]]] <- appended_info_species$CSize[[i]]
    A$logCsize[[appended_info_species$taxa[i]]] <- appended_info_species$logCsize[[i]]
    A$coords[[appended_info_species$taxa[i]]] <- as.array(appended_info_species$coords[[i]])
    A$species[[appended_info_species$taxa[i]]] <- appended_info_species$species[[i]]
  }
}
A$taxa <- names(A$CSize)

#######################################################################################################################################

# STEP 3: Perform pairwise comparisons 

# Pairwise comparisons test whether species differ in ontogenetic trajectories

# Iteratively generating Procrustes linear model for each species and save slope/intercept coefficients
TrajectoryList=list()
SlopesList=list()
InterceptsList=list()
ElevationsList=list()
for (j in 1:length(A$taxa)){
  x<-procD.lm(A$PCvalues[[j]]~A$logCsize[[j]],iter=999) # Using A$coords or A$PCvalues gives the same output for the pairwise comparisons, need to use A$PCvalues here to plot the trajectories
  TrajectoryList[[A$taxa[[j]]]]<-x
  SlopesList[[A$taxa[[j]]]]<-x$coefficients[2,]
  InterceptsList[[A$taxa[[j]]]]<-x$coefficients[1,]
  ElevationsList[[A$taxa[[j]]]]<- (x$coefficients[2,] * mean(log(A$CSize[[j]])) ) + x$coefficients[1,]
}


# Iteratively create two species pairs and perform Procrustes ANOVA and save output data from comparison
# Pairwise Procrustes ANOVAs of ontogenetic trajectories
PairwiseComparisons=list()
SlopeDifferences=list()
InterceptDifferences=list()
g=1
n=2
for (g in 1:length(A$taxa)){
  for (h in n:length(A$taxa)){
    if (g==length(A$taxa)){
      break
    }else
      
    SpeciesA <- A$taxa[g]
    SpeciesB <- A$taxa[h]
    toMatch <- c(SpeciesA, SpeciesB)
    filename <- paste(SpeciesA, SpeciesB, sep=" vs. ")
    
    two.species <- abind(A$coords[[SpeciesA]],A$coords[[SpeciesB]])
    two.species.csize <- c(A$CSize[[SpeciesA]],A$CSize[[SpeciesB]])
    two.species.species <- as.factor(c(paste(A$species[[SpeciesA]]),paste(A$species[[SpeciesB]])))
    two.species.slope.diff<-SlopesList[[g]]-SlopesList[[h]]
    names(two.species.slope.diff) <- paste("PC", c(1:164), " slope diff.", sep="")
    two.species.intercept.diff<-InterceptsList[[g]]-InterceptsList[[h]]
    names(two.species.intercept.diff) <- paste("PC", c(1:164), " int. diff.", sep="")
    
    Output <- procD.lm(two.species~two.species.csize, ~two.species.species, logsz = TRUE, iter = 9999)##Remember to add more iterations##
    # logsz above = TRUE which means this is done for logged CS data
    PairwiseComparisons[[filename]]<-Output
    SlopeDifferences[[filename]]<-two.species.slope.diff
    InterceptDifferences[[filename]]<-two.species.intercept.diff
    
  }
  n=n+1
}


# Create matrix and fill with values from pairwise comparisons
# Fills using PairwiseComparisons object using $aov.table - pairwise Procrustes ANOVA of ontogenetic trajectories
Pvalues = matrix (nrow = length(names(PairwiseComparisons)),
                  ncol = length(c(paste("AOV", names(PairwiseComparisons[[g]]$aov.table[1,]), sep = "_"), names(SlopeDifferences[[g]][1]), names(SlopeDifferences[[g]][2]), names(InterceptDifferences[[g]][1]), names(InterceptDifferences[[g]][2]))),
                  dimnames = list(c(names(PairwiseComparisons)), c(paste("AOV", names(PairwiseComparisons[[g]]$aov.table[2,]), sep = "_"), names(SlopeDifferences[[g]][1]), names(SlopeDifferences[[g]][2]), names(InterceptDifferences[[g]][1]), names(InterceptDifferences[[g]][2])))
)
for (i in 1:length(names(PairwiseComparisons))){
  Pvalues[i,] <- rbind(as.numeric(paste(c(PairwiseComparisons[[i]]$aov.table[1,], SlopeDifferences[[i]][1], SlopeDifferences[[i]][2], InterceptDifferences[[i]][1], InterceptDifferences[[i]][2]))))
}


# Bonferroni correction to account for pairwise comparisons
# Create an empty matrix to fill later
CorrectedPvalues <- matrix(nrow = length(names(PairwiseComparisons)),
                           ncol = length(c("corrected AOV p.value")),
                           dimnames = list(c(names(PairwiseComparisons)), 
                                           c("corrected_AOV_p.value")))

# Fill the matrix with the Bonferroni corrected p-values
CorrectedPvalues[,1] <- p.adjust(Pvalues[,"AOV_Pr(>F)"], method = "bonferroni", 
                                 n = length(PairwiseComparisons))


#######################################################################################################################################

# STEP 4: Format the pairwise comparison results for plotting

# Convert pairwise results to a dataframe and add necessary columns
Pvalues <- as.data.frame(Pvalues)
CorrectedPvalues <- as.data.frame(CorrectedPvalues) %>% mutate(AOV_Rsq = Pvalues$AOV_Rsq)
CorrectedPvalues_sig <- CorrectedPvalues %>% mutate(sig_p = ifelse(corrected_AOV_p.value < .05, T, F),
                                                  p_if_sig = ifelse(sig_p, corrected_AOV_p.value, NA),
                                                  Rsq_if_sig = ifelse(sig_p, AOV_Rsq, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 3)))

# Add the pairwise comparison species names to seperate columns 
# This enables the plotting of a heatmap for pairwise differences
rownames <- as.character(rownames(CorrectedPvalues_sig))
rn <- do.call(rbind, strsplit(rownames, 'vs.'))
CorrectedPvalues_sig <- CorrectedPvalues_sig %>% mutate(Var1 = rn[,1], Var2 = rn[,2])
CorrectedPvalues_sig$Rsq_if_sig <- round(CorrectedPvalues_sig$Rsq_if_sig, digit =2)

# Create a palette for plotting the heatmap plot
mypalette_seq <- brewer.pal(9,"Oranges")
image(1:9,1, as.matrix(1:9), col = mypalette_seq,xlab="Oranges (sequential)",
      ylab = "", yaxt = "n")


#######################################################################################################################################

# STEP 5: Plot the results

# Plot a heatmap for the pairwise comparison
# Darker oranges indicate greater significance
Allometry_pairwise_heatmap <- ggplot(data = CorrectedPvalues_sig, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq[9], high = mypalette_seq[2], mid = mypalette_seq[5], 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", 
                       na.value =  mypalette_seq[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Ontogenetic allometry pairwise differences")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(face = "italic", angle = 45, size = 14, vjust = 1, hjust = 1),
        axis.text.y =  element_text(face = "italic", size = 14, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.5,
                               title.position = "top", title.hjust = 0.5))
Allometry_pairwise_heatmap


#####

# Plot a histogram linking the altricial-precocial spectrum with the number of ontogenetic trajectory pairwise differences

# Read in the data
# This data contains info for the altricial/precocial stage and number of sig differences from ontogenetic allometric trajectories
spectrum_data <- read.csv("Data/Precocial_altricial_proportion_of_sig_diffs_allometry.csv")
# Fix the species order to be the same as the .csv file
species <- factor(spectrum_data$Species,levels=unique(spectrum_data$Species))
# Create a colour palette
col.dev <- c("bisque", "sandybrown", "tan4", "lightskyblue3", "dodgerblue4")

# Plot the histogram
# Histogram of altricial-precocial spectrum and relationship to the number of significant differences from allometric trajectories
a_p_plot <- ggplot(spectrum_data, aes(colour = Group, fill = Group))+ 
  geom_bar(stat = "identity", aes(x=species, y = Proportion_of_sig_diffs_from_allometry))+
  theme_classic()+
  ylim(0,20)+
  xlab("Altricial-precocial spectrum")+ 
  ylab("Proportion of allometric significances")+
  ggtitle("Altricial-precocial spectrum and allometry relationship")+
  scale_color_manual(name = "Altricial-precocial spectrum", 
                     labels=c("Group_1" = "Super altricial", "Group_2" = "Altricial",
                              "Group_3" = "Semi-altricial", "Group_4" = "Semi-precocial", 
                              "Group_5" = "Precocial"),
                     values = col.dev, aesthetics = c("color","fill"))+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), 
        axis.text.x =  element_text(face = "italic", angle = 45, size = 10, vjust = 1, hjust = 1), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(vjust = -0.5, size = 12),  axis.title.y = element_text(hjust = 0.5, size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 13))
a_p_plot


#####

# Plot the trajectory output
# log(CS) vs PC1 with regression line plotted for each species

# The PC axes you want to plot
PC=1
# Determine x and y axis limits
Xlim<-1.1*c(min(log(gdf$CSize)),max(log(gdf$CSize))) # The values for the x-axis are always the same for each PC
Ylim<-1.1*c(min(PC_scores_PC1_flipped[,PC]),max(PC_scores_PC1_flipped[,PC]))
# Initiate the plot
plot(0, 0, type = "n",
     xlim = c(2,max(Xlim)),
     ylim = c(-0.3,max(Ylim)), # Adjust the value here for the most negative Ylim value
     xlab = "log(Centroid Size)",
     ylab = xlab(paste0("PC1 (", signif((PCA_summary[2,1]*100),4), "% of total variance)")), # Adjust to read the correct PC number on axis legend
     axes = FALSE,
     frame.plot = FALSE,
     asp=F)
# Plot the axes and zero-line
axis(1, c(2,3,4,5,6,7,8), pos=-0.3)
axis(2, c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), pos=2)
clip(-0.3,max(Xlim),-0.3,max(Ylim))
abline(h=0, lty=3)
# Plot the points for each specimen
points(log(gdf$CSize), PC_scores_PC1_flipped[,PC], pch=21, bg=alpha(col.species, 0.75), asp=F)
# Plot the ontogenetic trajectory regression lines
for (i in 1:length(TrajectoryList)){
  Line_coefficients <- TrajectoryList[[i]]$coefficients[,PC]
  Line_color <- mammal.PlottingValues$legendcolor[names(TrajectoryList)[i]]
  CS <- log(appended_info_species$CSize[[match(names(TrajectoryList)[i], appended_info_species$taxa)]])
  clip(min(CS),max(CS),-0.3,0.3)
  abline(Line_coefficients, col=Line_color, lwd=2)
}


### END ###
