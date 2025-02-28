
# Clear the entire working environment

rm(list = ls()) # clean all loaded data in the working environment
dev.off() # clears the plots
cat("\014") # cleans the console


# Load the Libraries
require(here) 
library(pvclust) # provides p-values for hierarchical clustering based on multiscale bootstrap resampling
library(vegan)
library(phytools)
library(maps)
library(stats)
library(cluster)
library(tidyverse)
require(sf)
require(sp)
require(fasterize)
require(rasterize)
require(rgeos)
require(dismo)

require(dplyr)
require(ggplot2)
require(stringr) # use function to glue string

library(stars)
library(rgdal)
library(raster)
library(terra)

library(factoextra) # necessary for conduction gap analysis
library(cluster)
library(vegan) # for performing hierachical clustering

library(pvclust) 
require(stringr) # use function to glue string for setting CRS
require(rgeos)
require(dismo)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### Apply Bootstrap HCA to PLANT SURVEY DATA
### Assessment of Stability Within the Cluster Groups
## Making Reference to the Optimal Number of K (clusters) from Gap Statistics

setwd("/Users/fnumbisi/Documents/RBGE_MajorFloras_CMEP/WADI_ASHAR_Project/WadiAshar_Sept2022Data_Clusters/GailClusterAnalysisDescription/AsharDAFORscale/DAFORspeciesBootstrapClusters")

# setwd("/Users/fnumbisi/Documents/RBGE_MajorFloras_CMEP/ALULA_FLORA_Project/AlUla_VegMapping_Data/AlUla_ClusterAnalysis_PlantAssemblages/Combined_versions")

# load the fpc package
#install.packages("fpc")
library(fpc) 

# read csv file for species list
# ------------------------------


AsharSpeciesDAFOR <- data.frame(read.csv(file.choose("Select CSV SpeciesList for 120 Survey points"), header=TRUE, sep = "," )) 

View(AsharSpeciesDAFOR[, 145:185])

AsharDAFORdata <- AsharSpeciesDAFOR[,1:183] # Subset the DAFOR Data to include only the species columns

rownames(AsharDAFORdata) <- AsharDAFORdata[,1] # rename row using the first column (PlotID)

AsharDAFORdata <- AsharDAFORdata[-1] # Remove the plot ID Columns

# Extract the needed data for Agglomerative Hierarchical Clustering (HCA)

#   Use all the columns except the first 
#   
# Scale the data columns to be zero mean and unit variance.
# The output of scale() is a matricx.

MetricsVar_scaled <- scale(AsharDAFORdata)   

#remove rows with missing values
MetricsVar_scaled <- na.omit(MetricsVar_scaled)

rownames(MetricsVar_scaled) <- AsharSpeciesDAFOR$PlotID # rename rows of scaled data



#### HIERARCHICAL CLUSTERING.
### Create dissimilary matrix.
# Use vegdist function to create dissimilarity) matrix using bray curtis method.
Ashar_bray <- vegdist(AsharDAFORdata, method = "bray", binary = FALSE)
# Needs to be binary = TRUE for presence/absence dataset.

### Cluster analysis. 
# Create a cluster result object from the dissimilarity matrix created above, using the hclust command and select a linkage method.
# UWPMGA ("average") seems to be appropriate for our data but doesn't produce clear groups.
Ashar_bray_wardD <- hclust(Ashar_bray, method = "ward.D")

# Plot the objects as dendrograms.
plot(Ashar_bray_wardD, cex = 0.6)
str(Ashar_bray_wardD)


##----------GAP ANALYSIS
## Evaluate Optimum Number of Clusters For the scaled data

set.seed(20)
gap_stat_MetricsVar_scaled <- clusGap(MetricsVar_scaled, FUN = hcut, nstart = 25, K.max = 10, B = 50)

fviz_gap_stat(gap_stat_MetricsVar_scaled) #plot of clusters vs. gap statistic





#   Create the distance matrix.
# MetricsVar_disMat <- dist(MetricsVar_scaled, method="euclidean") # Use Euclidean Distance

# Get the HCA based on matrix of scaled data, Euclidean Distance, and Ward's Method
# Before running clusterboot() we’ll cluster the data using a hierarchical clustering 
# algorithm (Ward’s method):
#   Run the Hierarchical Clustering on the distance matrix. 
# MetricsVar_EuClust <- hclust(MetricsVar_disMat, method = "ward.D2")   


# The Gap Analysis suggests six clusters as shown on the gap statistics figure
# You can draw the rectangles on the dendrogram using the command 

#   Plot the dendrogram.
#plot(MetricsVar_EuClust)    
#rect.hclust(MetricsVar_EuClust, k=6)




#----------------
# B) Run the Bootstrap Clustering on the Topographic Variables Matrix 
# Let’s run clusterboot() on the protein data, using hierarchical clustering with five clusters.
# load the fpc package
# install.packages("fpc")
library(fpc)   

# set the desired number of clusters                               
kbest.p <- 9   # number of clusters obtained from gap analysis (statistics) curve    

#   Run clusterboot() with hclust 
#   ('clustermethod=hclustCBI') using Ward's method 
#   ('method="ward"') and kbest.p clusters 
#   ('k=kbest.p'). Return the results in an object 
#   called cboot.hclust.
# NB: By default clusterboot() runs 100 bootstrap iterations.

# NB: clusterboot() has a random component, so the exact stability values 
## and number of times each cluster is dissolved will vary from run to run.
# We reset the bootstrap iteration to 1000 r


set.seed(120) # set seed to similar random values in resampling
DaforVars_cboot.hclust <- clusterboot(Ashar_bray, B=1000, clustermethod = hclustCBI,
                                     method="ward.D2", k = kbest.p)

str(DaforVars_cboot.hclust)


#   The results of the clustering are in 
#   cboot.hclust$result. The output of the hclust() 
#   function is in cboot.hclust$result$result. 
#
#   cboot.hclust$result$partition returns a 
#   vector of clusterlabels. 

BoostrapClustergroups <- DaforVars_cboot.hclust$result$partition  # Extract the cluster groups or partitions

table(BoostrapClustergroups)
DaforDendrogramData <- DaforVars_cboot.hclust$result$result

k234 <- cutree(DaforDendrogramData, k = 9)

str(DaforDendrogramData)
plot(DaforDendrogramData)
boxcols <- c("sandybrown", "green", "darkred", "blue", "yellow", "purple", "red", "cyan",
          "orange")

cluster_type <- factor(DaforVars_cboot.hclust$partition)

rect.hclust(DaforDendrogramData, k=9, border = boxcols, cluster = k234)


library(dendextend)
dend1 <- color_branches(dend, k = 9, groupLabels = TRUE)
plot(dend1)

legend("topright", legend = levels(k234), fill = boxcols)



# install.packages("dendextend")
install.packages("circlize")
library(dendextend)
library(circlize)

# Create the dendrogram:
dend <- as.dendrogram(DaforDendrogramData)

# extra: showing the various clusters cuts 
k234 <- cutree(dend, k = 9)

# color labels by car company:
labels_colors(dend) <- DaforDendrogramData$labels[order.dendrogram(dend)]
# color branches based on cutting the tree into 9 clusters:
dend <- color_branches(dend, k = 9)


#plot(dend, horiz = TRUE)
plot(dend)
colored_bars(cbind(k234[,3:1], col_car_type), dend, rowLabels = c(paste0("k = ", 4:2), "Car Type"), horiz = TRUE)
legend("topright", legend = levels(car_type), fill = cols_4)




# --------------------------------------------------------------------------
# -- Run Bootstrap Clustering Results --------------------------------------

#   Create a function to print out the 
#   survey plots in each cluster, along with the values 
#   
#   
print_clusters <- function(dataVars,labels, k) {             
  for(i in 1:k) {
    print(paste("cluster", i))
    print(dataVars[labels==i,c(colnames(dataVars))])
  }
}


##---------------------------------------------------------------------------
# FUNCTION To create the dataframe for each cluster and append to list of dataframes

data_clusters <- function(DataFrame, labels, k) {             
  bootcluster_grp <- list()
  for(i in 1:k) {
    bootclust_id <- paste("cluster", i)
    bootclust_data <- DataFrame[labels==i,c(colnames(DataFrame))]
    bootclust_data <- cbind(bootclust_id, bootclust_data)
    bootcluster_grp[[bootclust_id]] <-bootclust_data
  }
  return(bootcluster_grp)
}


# --- Use the function ("print_clusters") to print the data for the points in each cluster

print_clusters(AsharDAFORdata, BoostrapClustergroups, kbest.p)



table(BoostrapClustergroups)

DaforVars_cboot.hclust$result # Get the results of the bootstrap clustering
str(DaforVars_cboot.hclust)


# --- Use the function ("data_clusters") to extract the Species Data for the points in each cluster
# Here, we use the cluster groups from the HCA with 1000 Bootstrap Resampling
bootClustersData2 <- data_clusters(AsharDAFORdata, BoostrapClustergroups, kbest.p)


#Write list of data.frames to separate CSV files with lapply() function
sapply(names(bootClustersData2), 
       function (x) write.table(bootClustersData2[[x]], file=paste(x, "BootstrapClusterSummaryDAFORSpecies.csv"), append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE) )  




# ------------------- BOOTRESULT --------------------------------------
# GET Matrix of Jaccard similarities for bootmethod="boot". 
# Rows correspond to clusters in the original data set. 
# Columns correspond to bootstrap runs.

JaccardMat <- DaforVars_cboot.hclust$bootresult # Get matrix of Jaccard similarities  

## -------------------------------------
### Evaluate Stability within Clusters

# The vector of cluster stabilities. 
# Values close to 1 indicate stable clusters

DaforVars_cboot.hclust$bootmean  # 	Get clusterwise means of the bootresult.                                 


# The count of how many times each cluster was 
# dissolved. By default clusterboot() runs 100 
# bootstrap iterations. We reset the iterations to 1000 
# Clusters that are dissolved often are unstable. 

DaforVars_cboot.hclust$bootbrd  # Get clusterwise number of times a cluster has been dissolved.                                  


##----CONFIDENCE INTERVAL OF JACCARD 
# Compute a 95% confidence interval for the mean Jaccard Similarity of the bootstrap clusters.


# append mean Jaccard similarity
JaccardMat_data <- cbind(JaccardMat, JSmean = DaforVars_cboot.hclust$bootmean) 


# Calculating standard deviation of each row
#  use apply and transform functions
JSstd <- as.data.frame(apply(JaccardMat,1, sd, na.rm = TRUE))
colnames(JSstd) <- "JSstd"

JaccardMat_data <- as.data.frame(append(as.data.frame(JaccardMat_data), JSstd))


# Calculate the margin of error for 95% CI
n=1000 # Set the number of iterations
JaccardMat_data$JSmargin <- qt(0.975,df=n-1)*JaccardMat_data$JSstd/sqrt(n)


# Determine the lower and upper confidence interval boundaries.

JaccardMat_data$JSlowerCI <- JaccardMat_data$JSmean - JaccardMat_data$JSmargin # Lower bound of 95% CI
JaccardMat_data$JSupperCI <- JaccardMat_data$JSmean + JaccardMat_data$JSmargin # Lower bound of 95% CI



#---------------------------------------------------------------------
# CLUSTER STABILITY ASSESSMENT
# Boxplot of Mean Stability index and 95% CI

DaforBootClusters <- c("Cluster1","Cluster2","Cluster3","Cluster4",
                  "Cluster5", "Cluster6", "Cluster7", "Cluster8", "Cluster9")

# Append column "BootClusters" to Jaccard Similarities Dataframe
JaccardMat_data2 <- cbind(JaccardMat_data, DaforBootClusters = as.factor(DaforBootClusters)) 

write.table(JaccardMat_data2, file= "HCABoostrapJaccardSimilaritesDAFORspecies.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)

#install.packages("plotrix")
library(plotrix)
library(ggplot2)

# ggplot2 plot with confidence intervals

ggplot(JaccardMat_data2, aes(DaforBootClusters, JSmean)) +
  geom_point() +
  geom_errorbar(aes(ymin = JSlowerCI,
                    ymax = JSupperCI)) +
  xlab("Clusters \n (Bootstrap resampling, 1000 iterations)") +
  ylab("Jaccard Similarity index")




##### APPEND BOOTSTRAP CLUSTER GROUPS TO SPECIES DATAFRAMES
## --------------------------
### Append the Bootstrap Cluster Groups to Species List and Survey Points Dataframe

AsharSpeciesBin <- AsharSpeciesDAFOR # Duplicate species list


# Convert Species list from DAFOR to Binary
AsharSpeciesBin[2:183][AsharSpeciesBin[2:183] > 0] <- 1 # replace all DAFOR codes > 0 to 1


AsharSpeciesPts <- read.csv(file.choose("Select CSV Surved Descriptn for 120 Survey points"), header=TRUE, sep = "," ) 


BINARYSpeciesClusters <- cbind(AsharSpeciesBin, DaforBootCluster = BoostrapClustergroups)

DAFORSurveyPtsClusters <- cbind(AsharSpeciesPts, DaforBootCluster = BoostrapClustergroups)

DAFORSpeciesBootstrapClusters <- cbind(AsharSpeciesDAFOR, DaforBootCluster = BoostrapClustergroups)

#DAFORSpeciesBootstrapClusters <- cbind(DAFORSpeciesBootstrapClusters, DAFORCluster = AsharSurveyPtsDesc2$DAFORscale_group)


write.table(DAFORSpeciesBootstrapClusters, file= "DAFORSpeciesBootstrapClustersUpdate.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


write.table(BINARYSpeciesClusters, file= "BINARYSpeciesBootstrapClustersUpdate.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
write.table(DAFORSurveyPtsClusters, file= "DAFORSurveyPtsBootstrapClusters.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
##-----TALLY SPECIES IN EACH SURVEY POINT
# Get row-wise sum of Species column and append (column bind) to dataframe

BINARYSpeciesClustersSUM <- data.frame(cbind(BINARYSpeciesClusters, data.frame(sumSpecies = apply(BINARYSpeciesClusters[2:183], 1, sum))))

write.table(BINARYSpeciesClustersSUM, file= "BINARYSpeciesClustersSUM.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)







#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#----- SPLIT THE DATAFRAMES CLUSTERWISE
## Split the Species Dataframe by the number of cluster in the column "BootCluster"

clustList1 <- split(BINARYSpeciesClustersSUM, BINARYSpeciesClustersSUM$DaforBootCluster)

# Specify names of new data frame for each CLUSTER
split_names <- c("DafCluster1", "DafCluster2","DafCluster3","DafCluster4",         
                 "DafCluster5", "DafCluster6", "DafCluster7", "DafCluster8", "DafCluster9")

for (i in 1:length(clustList1)) {        # Run for-loop
  assign(split_names[i], clustList1[[i]])
}


# Create a cross tabulation that displays the frequency distribution of TopoCluster per DAFORclusters
# Convert the output back into a data frame.
# 
#
#CLusterTab <- table(BINSurveyPtsClusters$BootCluster, BINSurveyPtsClusters$DAFORscale_group)
#CLusterCompared <- as.data.frame.matrix(CLusterTab)


#CLusterCompared2 <- as.data.frame.matrix(table(BINSurveyPtsClusters$DAFORscale_group, BINSurveyPtsClusters$BootCluster))






#-----------------------------------------------------------------
#------- FUNCTION FOR EXTRACTING SPECIES FOR EACH CLUSTER

library(caret)
library(dplyr)

SpeciesInCluster <- function(clusterList){
  
  require(caret) # library needed to use the "select" function
  require(dplyr)
  
  colStart <- as.integer(readline(prompt = "Enter index number of first species column in dataframe: "))
  colEnd <- as.integer(readline(prompt = "Enter index number of last species column in dataframe: "))
  
  clusterName <- deparse(substitute(clusterList)) # extract the name of dataframe.
  
  # Extract columns with non-zero column sum
  #cluster_1a <- clusterList[colStart:colEnd] %>% select(which(!colSums(clusterList[colStart:colEnd], na.rm=TRUE) %in% 0))
  cluster_1a <- clusterList[colStart:colEnd] %>% select_if(~ !is.numeric(.) || sum(.) != 0)
  
  # Compute sum of Columns
  cluster1Sum <- data.frame(colname = names(cluster_1a), sumPlotOccurrence=colSums(cluster_1a))
  
  cluster1Sum <- data.frame(cluster1Sum[order(-cluster1Sum$sumPlotOccurrence), ])
  rownames(cluster1Sum) <- NULL
  colnames(cluster1Sum)[1] <- "Species"
  
  FinclustName <- paste(clusterName, "Species", sep = "") # set name for cluster i dataframe
  
  FinclustName <- cluster1Sum 
  
  write.table(FinclustName, file= paste(clusterName, "_Bootstrap_Species.csv", sep = ""), append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
  
  return(FinclustName)
}


# Use the function "SpeciesInCluster" to extract species in plots within each cluster
# and tally the number of plots in which they occur.

SpeciesInCluster(DafCluster9)



##------------------------------------------------------------
##------------------------------------------------------------

#  INDICATOR SPECIES ANALYSIS
#install.packages("labdsv")

library(labdsv)
library(indicspecies)

#### INDICATOR SPECIES.
### Indicatorspecies analysis.
# Use multipart function on species x sites data matrix (AsharSpeciesBin[2:183]) and groups vector (BoostrapClustergroups)
# to determine lists of species that are associated to grouped sites/samples).


Ashar_BootIndicatorSpeceis <- multipatt(AsharSpeciesDAFOR[2:183], BoostrapClustergroups, 
                                        control = how(nperm=999)) 

# Display the indicator value components ‘A’- the specificity or positive predictive value and 
#  ‘B’ - the fidelity or sensitivity of the species as indicator of the target site group

# "A" is the probability that a site is in the target group given that the species has been found. 
# It is the mean abundance of species in the target group divided by the sum of the 
# mean abundance values over all groups.

# "B" is the probability that you will find the species in a site that is in the target site group. 
# It is the relative frequency of occurrence of the species in the target group.
# High B value means that if you do not find the species at a given site, 
# the probability that the site belongs to the target group is low.

summary(Ashar_BootIndicatorSpeceis, indvalcomp=TRUE)


# Get a more comprehensive list of species and group associations, 
# including ones that were not significant, by changing the significance level:
summary(Ashar_BootIndicatorSpeceis, indvalcomp = TRUE, alpha = 1)

# For some species, the indicator value is the highest when the groups are combined. 
# Use the object "$sign" to find these species.
Ashar_BootIndicatorSpeceis$sign




