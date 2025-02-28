
# Clear the entire working environment

rm(list = ls()) # clean all loaded data in the working environment
dev.off() # clears the plots
cat("\014") # cleans the console


#install.packages('devtools') #assuming it is not already installed
#library(devtools)


version
packageStatus()
## ------------------------------------------------------

#setwd('/Users/fnumbisi/Documents/RBGE_MajorFloras_CMEP/ALULA_FLORA_Project/AlUla_AOI_Data')

setwd('/Users/fnumbisi/Documents/RBGE_MajorFloras_CMEP/ALULA_FLORA_Project/')

## ------------------------------------------------------
# Load the necessary libraries

#knitr::opts_chunk$set(echo = TRUE) # use when running rmarkdown

# Load the Libraries


library(here) 
library(rgdal) 
library(tidyverse)
library(raster)
library(sf)
library(sp)
library(fasterize)
library(rasterize)
library(rgeos)
library(dismo)

library(dplyr)
library(ggplot2)
library(stringr) # use function to glue string

#install.packages("snowfall")

## ------------------------------------------------------

list.of.packages <- c("plyr", "parallel", "ranger", "raster",
                      "rgdal", "rgrass7", "snowfall", "lidR", "knitr", "tmap")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages, dependencies = TRUE)}


## ------------------------------------------------------
library(rgdal)
library(raster)
library(terra)
#install.packages("plotKML")
#library(plotKML)
library(sp)

#install.packages("devtools")
#library(devtools)
library(snowfall)
library(parallel)
library(doParallel)


#install.packages("GSIF")
#library(GSIF)

GDALinfo("./AlUla_AOI_Data/2023AlUlaIPVI2stdIQR.tif")

# Create Stack of VIs STD for whole of AlUla

AlUlaVIsStd <- raster::stack(paste0("./AlUla_AOI_Data/2023AlUla", c("NDVI","IPVI2","SAVI","SARVI2","GSAVI","GRVI","MSAVI2"), "stdIQR.tif"))

plot(AlUlaVIsStd)

canProcessInMemory(AlUlaVIsStd, verbose = TRUE) #test if raster object can be loaded into memorr for precessing

# Caculate memory requirements assuming 8 bytes per cell value
n_values <- ncell(AlUlaVIsStd) * nlayers(AlUlaVIsStd) # get number of values in rasterstack
mem_estimate <- (8 * n_values) / 2^20 # get memory estimate used by rasterstack (in mb)


## Use Parallelisation to estimate actual memory usage by rasterstack

# Start a cluster with 1 less the number of cores in 

clust <- makeCluster(parallel::detectCores() - 1)

registerDoParallel(clust) # register the parallel backend cluster

#system.time({mem_actual <- as.integer(object.size(readAll(AlUlaVIsStd))) / 2^20 # get actual memory usage by the rasterstack
#})

# stop cluster
#stopCluster(clust) # Vector Memory Exhausted (limi reached)


#mem_actual <- as.integer(object.size(readAll(AlUlaVIsStd))) / 2^20 # get actual memory usage by the rasterstack

print(c('Talal Values: ',n_values, 'Est. Memory Used: ',mem_estimate))# 'Actual Memory Used: ', mem_actual))


#-----------------------------------------------------------------------------------
# Will need to process the rasterstack in blocks of row or by parallelisation tiles, 
# each small enough to store in memory

# temporary output file
# rstk_out <- tempfile(fileext = ".tif")

# Identify the number of blocks (or tiles) to process the rasterstack
rstk_blocks <- blockSize(AlUlaVIsStd)
print(rstk_blocks)


## --------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

stack_list <- paste0("./AlUla_AOI_Data/2023AlUla", c("NDVI","IPVI2","SAVI","SARVI2","GSAVI","GRVI","MSAVI2"), "stdIQR.tif")

r <- raster::raster("./AlUla_AOI_Data/2023AlUlaIPVI2stdIQR.tif")

AlUlaVIsStd2 <- terra::rast(stack_list) # read and stack large raster using the terra package
class(AlUlaVIsStd2)

plot(AlUlaVIsStd2)
plot(AlUlaVIsStd2[[1]]) #GRVI
plot(AlUlaVIsStd2[[3]]) #IPVI
plot(AlUlaVIsStd2[[5]]) #NDVI
plot(AlUlaVIsStd2[[7]]) #SAVI


## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
AlUlaSites <- rgdal::readOGR(dsn="./AlUla_AOI_Data/AlUla_Polygons_04Oct2023_UTM_GS/", 
                             layer="AlUla_Polygons_04Oct2023_UTM") 

plot(AlUlaSites)

AlUlaSites_Data <- as.data.frame(AlUlaSites)


AlUlaSites_Data_sub <- AlUlaSites_Data[,1:8]
AlUlaSites_Data_sub
  


##------------------------
# Efficient spatial overlay and raster data extraction using terra package

# Check if files for identical CRS
identicalCRS(AlUlaVIsStd, AlUlaSites)

# transform CRS of reference sites CRS to R that of raster

AlUlaSites_LatLong <- sp::spTransform(AlUlaSites, CRSobj = "+proj=longlat +datum=WGS84 +no_defs")

proj4string(AlUlaVIsStd) # get the project of raster data


proj4string(AlUlaSites_LatLong)



crs(AlUlaSites_LatLong) <- proj4string(AlUlaVIsStd)

identicalCRS(AlUlaVIsStd, AlUlaSites_LatLong)



# Extract Raster values at point
AlUlaVIstdData <- terra::extract(rast(AlUlaVIsStd), # Use image stack based on raster package
                                 vect(AlUlaSites, crs="+proj=longlat +datum=WGS84 +no_defs"))






## ----------------------------------------------------------------------------
# Use Parallelisation to stack and tile VIs Images for AlUla







##-------------------------------------------------------

AlUla <- readGDAL("./AlUla_AOI_Data/2023AlUlaIPVI2stdIQR.tif")

#AlUla.tile.lst <- getSpatialTiles(obj, block.x=1000)


## raster files via rgdal:
#library(rgdal)
#fn = system.file("pictures/SP27GTIF.TIF", 
#                 package = "rgdal")
obj <- GDALinfo(AlUla)

obj2 <- GDALinfo("./AlUla_AOI_Data/2023AlUlaIPVI2stdIQR.tif")
obj2

obj3 <- GDALinfo("./AlUla_AOI_Data/AlUla_30msrtm_dem.tif")

# Read raster using the terra library
obj3 <- terra::rast("./AlUla_AOI_Data/AlUla_30msrtm_dem.tif")
plot(obj3b)

obj4 <- terra::rast("./AlUla_AOI_Data/AlUla_30mSlopeDegree.tif")
plot(obj4)



#AlUla.tile.lst <- tile(obj2, block.x=1000)
#AlUla.tile.lst <- tile(obj2, block.x=5400)





#AlUlaoffset <- c(AlUla.tile.lst$offset.y[1], AlUla.tile.lst$offset.x[1])

#AlUlaregion.dim <- c(AlUla.tile.lst$region.dim.y[1], 
#                AlUla.tile.lst$region.dim.x[1])


## read the first tile:
#AlUla_GTIF_T1 <- readGDAL(AlUla, offset=AlUlaoffset, 
#                        region.dim=AlUlaregion.dim)


#str(AlUla_GTIF_T1)



##-------------------------------------------------------
##-------------------------------------------------------
# Tiling Raster using the meteo and raster packages

#install.packages("meteo")
library(meteo)
library(raster)
library(CAST)
library(sf)
library(ranger)
library(plyr)
library(parallel)
library(doParallel)





#Identify your polygon boundary
#"./AlUla_AOI_Data/AlUla_BoundReserves/AlUla_ProjectBound.shp" #directory where the .shp is 
ppath <- "./AlUla_AOI_Data/AlUla_BoundReserves/" #directory where the .shp is 
pname <- "AlUla_ProjectBound" #filename of the actual .shp without the extension.

#Read your polygon boundary
ppoly <- readOGR(ppath,pname) #if this doesn't work theyn maybe you do not have the package loaded? 
plot(ppoly)
class(ppoly) # SpatialPolygonDataFrame

#Identify rectangular spatial extent
e <- extent(ppoly)
plot(e)

# Convert SpatRaster to a SpatVector of Polygons
terra::as.

# Convert extent of Tile into polygons using the rast within the Terra package
Ext_polyg1 <- terra::as.polygons(ext(tiles_spdf[[15]])) 
class(Ext_polyg1)

crs(Ext_polyg1) <- "epsg:4326" # Set the CRS for the tile extent
crs(Ext_polyg1, describe=TRUE)$area
plot(Ext_polyg1)

# coerce the extent to a SpatialPolygons object using Sp library
Ext_polyg2 <- as(e, "SpatialPolygons")
plot(Ext_polyg2, add=TRUE)

plot(ppoly)
#plot(Ext_polyg2, add=TRUE)
plot(Ext_polyg1, add=TRUE)

# Plot All the Tile in 
plot(ppoly)


tiles_spdf


TilesPolyg <- list()

for (i in 1:length(tiles_spdf)) {
  
  Ext_polyg <- terra::as.polygons(ext(tiles_spdf[[i]])) 
  crs(Ext_polyg) <- "epsg:4326" # Set the CRS for the tile extent
  crs(Ext_polyg, describe=TRUE)$area
  #plot(ppoly)
  #plot(Ext_polyg, add=TRUE)
  
  TilesPolyg[i] <- Ext_polyg
}


TilesPolygMerge <- do.call(merge, TilesPolyg)

# Create Polygon from all tiles


Ext_polygALL <- terra::as.polygons(ext(tiles_spdf)) 
crs(Ext_polygALL) <- "epsg:4326" # Set the CRS for the tile extent
crs(Ext_polygAll, describe=TRUE)$area

# Combine SpatVectors
tiles_spdfAll <- c(tiles_spdf)
class(tiles_spdfAll)

for (i in 6:9) {
  
  plot(tiles_spdf[[i]])
  
}

plot(Ext_polyg1)
crs(Ext_polyg1)

plot(as(extent(tiles_spdf[[5]]), "SpatialPolygons"), add=TRUE)



#Projection type and datum, respectively
#p1 <- "+proj=UTM"
#p2 <- "+datum=WGS84"
#projmap <- paste(p1,p2)

#Identify the list of images you wish to perform the loop
files <- list.files(ipath, pattern=".tif$", full.names=FALSE)
stopifnot(length(files)>0)

#add output directory
outfiles <- pasteO(opath, files) #paste0 forgoes any separator between the two vectors "opath" and "files". if you need to add a slash somewhere then you need to change it to the paste function or change your vectors above.

#Change extension
extension(outfiles) <- "tif" #if you want to make it .tif 

for (f in 1:length(files)){
  r <- raster(files[f])
  rc <- crop(r, e)
  rm <- mask(rc, ppoly)
  rw <- writeRaster(rm, outfiles[f], overwrite=TRUE)
}





for (i in seq(length(tiles_spdf))) {
  
  stopifnot(length(tiles_spdf)>0) # interrupt loop if length of list is not >0
  
  imageTile <- tiles_spdf[[i]]
  
  img_extent <- extent(imageTile) # Get the extent of tile
  
  if (img_extent )
    
    writeRaster(imageTile, filename = ofile, format='GTiff',
                options='COMPRESS=LZW', datatype='FLT4S', overwrite=TRUE)
  
    
  # convert sf multi-polygon data frame to large SpatialPolygonsDataFrame
  ng2_sp <- as(ng2, Class = "Spatial") 
  
  str(ng2_sp) ###
  plot(ng2_sp) ###
  
  
  # Set the projection CRS for the species occurrence grid to that of AOI shapefile
  projection(ng2_sp) <- projection(studyAOI)
  
  # intersect from raster package: 
  
  AOIintersectGrid <- intersect(ng2_sp, studyAOI) # clip the grid to boundary of AOI (Soqotra Mainland)
  
  AOIintersectGrid$Cellid <- 1:nrow(AOIintersectGrid) # Create Id values for all cells in intersect boundary
  
  ng2_sp$Cellid <- 1:nrow(ng2_sp) # Create Cell id for all square grids of the AOI extent
  
  # Export clipped grid as shapefile to specified directory
  ## Create Directory to Save Results
  
  if (!dir.exists("AOIShapefileGrid")) dir.create("AOIShapefileGrid", recursive = TRUE)
  
  dir.create("AOIShapefileGrid") ###
  
  
  
}





rv <- raster(v, crs=CRS(projargs), xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx)

writeRaster(rv, filename = ofile, format='GTiff',
            options='COMPRESS=LZW', datatype='FLT4S', overwrite=TRUE)





# Identify Tiles where the training data falls in
AlUla.tile.Use <- sp::over(AlUlaSites, AlUla.tile.lst, returnList = TRUE)





t2<- tiling(r,tilesize = 5000,overlapping = 0) # produces 24 tiles for AlUla: each tile has 5000 cells
t2_tbl <- as.data.frame(t2) # create tiles table


par()
plot(t2[[8]], col=grey.colors(255)); plot(t2[[8]])


AlUla_6tiles_ipvi <- list.files(t)


plot(t[[50]],col=grey.colors(255))
plot(t[[50]],col='viridis')

# Save computed tiles as separate file in designated folder



list.files("./AlUla_VIstd_tiles/")




## -----------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Extract Values from Raster Data FOR RANDOM FOREST CLASSIFICATION


library(raster)

# Load a raster file
FilePath="C:\\Temp3\\LandSat_WGS83.tif"
raster1 <- raster::raster(FilePath)

# Load a CSV that contains colmns of X and Y coordinate values

pointCoordinates <- read.csv("C:\\Temp3\\LandSatCoordinates.csv")
coordinates(pointCoordinates) <- ~ X+ Y # Convert to a coordinates object

# Extract the values using the coordinates
rasValue <- raster::extract(raster1, pointCoordinates)

# Add the pixel values back into the data frame
pointCoordinates <- data.frame(pointCoordinates,rasValue)




## -----------------------------------------------------------------------------

# LOAD LIBRARIES
# getting started: load necessary packages
require(rgdal)
require(raster)
require(rasterVis)
require(lattice)


# LOAD AND REPROJECT DEM RASTER TO UTM CRS
#   Import data - first a 30m DEM from NASA SRTM data (http://doi.org/10.5067/MEaSUREs/SRTM/SRTMGL1N.003) in geotiff format (though raster() will read any format that rgdal::readGDAL() recognizes), and then a point shapefile (a standard ESRI .shp, though rgdal::readOGR() can import many other sorts of spatial data also). 

# commands: raster::raster(), rgdal::readOGR() 
AlUlaDEM <- terra::rast("./AlUla_AOI_Data/AlUla_30msrtm_dem.tif")
rm(AlUlaDEM)
AlUlaDEM2 <- raster::raster("./AlUla_AOI_Data/AlUla_30msrtm_dem.tif")  # read raster

# It will be easier to work with projected data, so we'll project this to UTM using the appropriate proj4 string (http://proj4.org/index.html) for the CRS (Coordinate Reference System) that we want.
# Get the UTM Projected CRS for Wadi Ashar (or AlUla)
#spatialref <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"

AlUlaDEMutm <- projectRaster(AlUlaDEM2, 
                            crs="+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

# Have a quick look at this to make sure nothing has gone terribly wrong - with the raster package loaded typing the name of a raster object will give you summary data about that object.
AlUlaDEMutm  
plot(AlUlaDEMutm)

# LOAD AND SUBSET REFERENCE FIELD DATA (SHAPEFILE)
# read .shp (note that to read a shapefile, "the data source name (dsn= argument) 
# is the folder (directory) where the shapefile is, and the layer is the name of 
# the shapefile (without the .shp extension)" (from the rgdal::readOGR documentation))

AlUlaSites <- rgdal::readOGR(dsn="./AlUla_AOI_Data/AlUla_Polygons_04Oct2023_UTM_GS/", layer="AlUla_Polygons_04Oct2023_UTM") 
plot(AlUlaSites)
str(AlUlaSites)
class(AlUlaSites)

AlUlaSites_sub <- intersect() AlUlaSites

# project polygons to UTM (ZONE 37 in KSA)
AlUlaSites_utm <- spTransform(AlUlaSites, 
                             "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 

plot(AlUlaSites_utm)



## EXTRACT DATAFRAME FROM SF OBJECT
# Set the st_geometry property to NULL.

library(sf)

# system.file(package = "stats") # Gets the root of package 'stats'
# system.file() # The root of the 'base' package

f <- "AlUla_AOI_Data/AlUla_Polygons_04Oct2023_UTM_GS/AlUla_Polygons_04Oct2023_UTM.shp"
f_read <- sf::st_read(f)

# Read the shapefile from the root of sf packages
#nc <-  sf::st_read(system.file("shape/n.shp",
#                               package="sf"), quiet = TRUE)


AlUlaSites_utm_df <- AlUlaSites_utm %>% st_drop_geometry()
class(AlUlaSites_utm_df)

# Convert sp object to dataframe
AlUlaSites_utm_df <- as.data.frame(AlUlaSites_utm_df)
class(AlUlaSites_utm_df)





# subset points to eliminate sites of uncertain date - i.e., select from 'sites' only 
# those rows in which the 'period' column is "EIA" or "GalRom".

#AlUlaSites_sub <- sites[sites$period == "EIA" | sites$period == "GalRom",]  


# drop unused levels (not strictly necessary but will avoid messiness when plotting data later)
#sites_sub$period <- factor(sites_sub$period) 

# project points to UTM (ZONE 37 in KSA)
#sites_sub_utm <- spTransform(sites_sub, 
 #                            "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") 


# Check the file (note that it is now a Spatial Points Data Frame, and 
# typing its name will give you an object summary).  Note that there is 
# a 'type' field that we won't work with here, but which could be incorporated into 
# this kind of analysis, e.g., by further sub-setting or by grouping data when boxplotting.

#sites_sub_utm  


# COMPUTE DEM DERIVATIVES (SLOPE AND ASPECT)
#We can use the DEM to calculate DEM derivatives.

# commands: raster::terrain() 

AlUla_slope <- raster::terrain(AlUlaDEMutm, opt = 'slope', unit = 'degrees')  #calculate slope

AlUla_aspect <- raster::terrain(AlUlaDEMutm, opt = 'aspect', unit = 'degrees') #calculate aspect

#Have a quick look at these to see that the results make sense - they are now raster objects just like areaDEM and can be examined the same way.
AlUla_slope

plot(AlUla_slope)

plot(AlUla_aspect)
AlUla_aspect


plot(AlUla_slope)
plot(AlUlaSites_utm, add = TRUE, col="blue")




## STACK RASTER DATA
# commands: raster::stack()
terrainstack <- stack(AlUlaDEMutm, 
                      AlUla_slope, 
                      AlUla_aspect)

terrainstack # have a quick look at resulting object, which shows the number of layers and the min/max values we expect

## EXTRACT RASTER DATA
# commands: raster::extract()
# extract the mean values w/in a 250m radius around each site for each terrain variable

# extract the mean values w/in each polygon

AlUlaSites_vals1 <- raster::extract(terrainstack,
                      AlUlaSites_utm,
                      #buffer = 250,
                      fun = mean,
                      sf = TRUE) 

AlUlaSites_vals2 <- raster::extract(terrainstack,
                                   AlUlaSites_utm,
                                   #buffer = 250,
                                   fun = mean,
                                   sp = TRUE) 

class(AlUlaSites_vals)

#AlUlaSites_vals[AlUlaSites_vals$]
parallel::detectCores()
test <- function(cpu=parallel::detectCores())
  #
# Convert values to dataframe by remove geometry

st_geometry(AlUlaSites_vals) <- NULL
class(AlUlaSites_vals)

old.packages()
update.packages()

AlUlaSites_utm_df







#-------------------------------------------------------------------------------
# Central Limit Theorem (CLT) application in R - I
#We will take sample size=30, 50 & 500 samples=9000
#Calculate the arithmetic mean and plot the mean of sample 9000 times

s30 <- c()
s50 <- c()
s500 <- c()

n =9000

# Generate 9000 replicates of different sample sizes, with replacment, from the population 
for ( i in 1:n){
  s30[i] = mean(sample(data$Wall.Thickness,30, replace = TRUE))
  s50[i] = mean(sample(data$Wall.Thickness,50, replace = TRUE))
  s500[i] = mean(sample(data$Wall.Thickness,500, replace = TRUE))
}
par(mfrow=c(1,3)) # set
hist(s30, col ="lightblue",main="Sample size=30",xlab ="wall thickness")
abline(v = mean(s30), col = "red")

hist(s50, col ="lightgreen", main="Sample size=50",xlab ="wall thickness")
abline(v = mean(s50), col = "red")

hist(s500, col ="orange",main="Sample size=500",xlab ="wall thickness")
abline(v = mean(s500), col = "red")

par(mfrow=c(1,1)) # Reset the plot window to 1row, 1 column

# generate 150 random number from normal distribution (i.e. mean of 0, sd of 1)
sampleNorm >- rnorm(150, mean=0, sd=1) 








#-------------------------------------------------------------------------------
## EXTRACT RASTER VALUES FROM VEGETATION INDICES
# ------------------------------------------------------------------------------
#### Extract raster information from the Raster Stack


referenceSites <- MergedRefPolygons # For polygons are reference data

referenceSitesData <- data.frame(referenceSites)
referenceSitesData %>% count(ClassType, sort = TRUE)


#referencePoints <- MergedRefPoints # For points are reference data

# As there are several packages with the extract() function, and possibly conflict if these loaded libraries
# we specify to use the function from the raster package using the "::" notation 

extrVIGLCMs <- raster::extract(AsharVI_VIGLCMraster, referenceSites, df=TRUE)


#extrVIs <- raster::extract(AsharVIstack, referencePoints, df=TRUE)



# Create a column with unique value for buffer polygons (for each survey point) - "PolygonID" 

referenceSites$PolygonID <- 1:nrow(referenceSites) # for reference site polygons

#referenceSites$PointID <- 1:nrow(referenceSites) # for reference site Points (Scenario 3)


# Merge the Extracted raster data and reference survey sites

extrVIGLCMs <- merge(extrVIGLCMs, referenceSites, by.x="ID", by.y="PolygonID") # Use in cas of polygons

#extrVIs <- merge(extrVIs, referenceSites, by.x="ID", by.y="PointID") # use in case of points (scenario 3)

write.table(extrVIGLCMs, file=paste(Wdir, "VIGLCMoutput/VIGLCMrasters/AsharVI_VIGLCM32Grey7_7window.csv", sep = ""),
            append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)


dim(extrVIGLCMs) # verify the size of extracted, merge data
View(extrVIGLCMs)
extrVIGLCMs %>% count(ClassType, sort = TRUE)


# head(extrVIGLCMs[,c(1:26,107)]) # For data with with compute 4 GLCM Textures for each VI
# tail(extrVIGLCMs[,c(1:26,107)])

head(extrVIGLCMs[,c(1:41, 122)]) # For data with with compute 4 GLCM Textures for each VI
tail(extrVIGLCMs[,c(1:41, 122)])

str(extrVIGLCMs[,c(1:41, 122)])

# View mean,median,25th and 75th quartiles,min,max
summary(extrVIGLCMs[,c(1:41, 122)])



# Convert Response Variable to factor with several classes
extrVIGLCMs$ClassType <- as.factor(extrVIGLCMs$ClassType)


### --------------------------------------------------------
# Split Data for Training: 
# Using the polygon buffer for 120 survey point, we had 827 observation points in extracted, merged data
#In order to speed things up, for this tutorial we will reduce the data. 
#Therefore, from each training polygon only 75% of the pixels will be used for model training.
#Therefore, from each polygon 75% of the pixels are randomly drawn.

# We us the createDataPartition() function to split training and test data
# set.seed(120)
# trainids <- createDataPartition(extrVIs$ID, list=FALSE, p=0.75)
# trainDat <- extrVIs[trainids,]
# nrow(trainDat)

# Use the entire data points for Scenario 1, 2, and 3 as RF does Bagging and OOB Resampling
trainDat <- extrVIGLCMs 
nrow(trainDat)

#--------------------------------------------------------
## Model training
### Predictors and response


# predictorsVIGLCM <- c("EVI", "GSAVI", "MSAVI2", "NDVI", "SAVI",
#                 "EVI_Var", "GSAVI_Con", "GSAVI_Ent", "GSAVI_Var",
#                 "MSAVI2_Con", "MSAVI2_Ent", "MSAVI2_Var",
#                 "NDVI_Con", "NDVI_Ent", "NDVI_Var",
#                 "SAVI_Con", "SAVI_Ent", "SAVI_Var")

# # All Correlation Textures have NAs, thus excluded from predictors
predictorsVI7GLCM_7Window32Grey <- c("EVI", "GSAVI", "MSAVI2", "NDVI", "SAVI",
                                     "EVI_Mea", "EVI_Var", "EVI_Hom", "EVI_Con", "EVI_Dis", "EVI_Ent",
                                     "GSAVI_Mea", "GSAVI_Var", "GSAVI_Hom", "GSAVI_Con", "GSAVI_Dis", "GSAVI_Ent",
                                     "MSAVI2_Mea", "MSAVI2_Var", "MSAVI2_Hom", "MSAVI2_Con", "MSAVI2_Dis", "MSAVI2_Ent",
                                     "NDVI_Mea", "NDVI_Var", "NDVI_Hom", "NDVI_Con", "NDVI_Dis", "NDVI_Ent",
                                     "SAVI_Mea", "SAVI_Var", "SAVI_Hom", "SAVI_Con", "SAVI_Dis", "SAVI_Ent") 

# # 2ns Predictor Groups VI, GLCM Textures for Mean and Homogeneity
predictorsVI7GLCM_7Window32Grey_2 <- c("EVI", "GSAVI", "MSAVI2", "NDVI", "SAVI",
                                       "EVI_Mea", "EVI_Hom", 
                                       "GSAVI_Mea", "GSAVI_Hom",
                                       "MSAVI2_Mea", "MSAVI2_Hom", 
                                       "NDVI_Mea", "NDVI_Hom", 
                                       "SAVI_Mea", "SAVI_Hom") 

# # 3rd Predictor Groups: GLCM Textures for Mean and Homogeneity
predictorsVI7GLCM_7Window32Grey_3 <- c("EVI_Mea", "EVI_Hom", 
                                       "GSAVI_Mea", "GSAVI_Hom",
                                       "MSAVI2_Mea", "MSAVI2_Hom", 
                                       "NDVI_Mea", "NDVI_Hom", 
                                       "SAVI_Mea", "SAVI_Hom") 


response <- "ClassType"

# --------------------------------------------------------
### Model training
# We then train a Random Forest model to lean how 
# the classes can be distinguished based on the predictors

# Caret's train() function can be used to train the model.
# Before model training, we can specify some control settings using trainControl. 
# For hyperparameter tuning (mtry) as well as for error assessment, 
# we use a spatial 7-fold cross-validation. Therefore, the training data are split into 7 folds 
# but data from the same polygon are always grouped so that they never occur in both 
#  training and testing. Also we make sure that each fold contains data from each land cover class. 
# For this, we use the CAST's CreateSpacetimeFolds() function for spatial cross-validation when
#  we specify the polygon ID and the class label ("ClassType").

set.seed(120)
indices <- CreateSpacetimeFolds(trainDat, spacevar = "ID", k=7, class="ClassType")
str(indices)
ctrl <- trainControl(method="cv", 
                     index = indices$index,
                     savePredictions = TRUE)


# train the model
set.seed(120)
model <- ffs(trainDat[,predictorsVI7GLCM_7Window32Grey_3],
             trainDat[,response],
             method="rf",
             metric="Kappa",
             trControl=ctrl,
             importance=TRUE,
             ntree=250)


# Verify the Model
print(model)
plot(varImp(model))



## -------------------------------------------------
### Model validation

# When we print the model (see above) we get a summary of the 
# prediction performance as the average Kappa and Accuracy of the three spatial folds. 
# Looking at all cross-validated predictions together we can get the "global" model performance.

# get all cross-validated predictions:
cvPredictions <- model$pred[model$pred$mtry==model$bestTune$mtry,]

# Calculate cross table: 
# Classificaton confusion Matrix: row = predictiong, columns = reference/0bserved classes
table(cvPredictions$pred,cvPredictions$obs)

# create confusion matrix using the consusionMatrix() function in Caret package
library(caret)
confusionMatrix(cvPredictions$pred, cvPredictions$obs,
                mode = "everything",
                positive="1")

confusTable <- data.frame(table(cvPredictions$pred,cvPredictions$obs))

write.table(confusTable, file= paste0(Wdir, "VIGLCMoutput/VIGLCMrasters/AsharVI_VIGLCM32Grey7_7window_ConfusionMatrix.csv"), 
            append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
### Model prediction
names(AsharVI_VIGLCMraster)
AsharVI_VIGLCMraster[["GSAVI_Hom"]]
# To perform the classification we can use the trained model and 
# apply it to each pixel of the raster stack using the predict function. Then we can 
# then create a map with meaningful colors of the predicted land cover using the tmap package.

VI7GLCM7Wind32GrePprediction <- predict(AsharVI_VIGLCMraster, model)

#cols <- c("sandybrown", "green", "darkred", "blue", "forestgreen", "sienna", "red")

# "black", "magenta", "maroon", "turquoise", "cyan4", "lightgreen", "maroon", "orange", "plum", "khaki"

cols <- c( "turquoise", "sandybrown", "magenta", 
           "yellow", "purple", "cyan4", "sienna", "green", "black", "red")

kmlColors2 <- c("#000201", "#72cea5", "#e69428", "#047c68", "#ec4f16", "#e8d22c", "#885bb5", 
                "#7d8b8f", "#3cec2c", "#91522d", "#d05e83")

# tm_raster(palette = cols, title = " Wadi Ashare \n Vegetation Assemblages")

tm_shape(deratify(VI7GLCM7Wind32GrePprediction)) +
  tm_raster(palette = kmlColors2, title = " Veg/LUC")+
  tm_scale_bar(bg.color="white", bg.alpha=0.75)+
  tm_layout(legend.bg.color = "white",
            legend.bg.alpha = 0.75)



## -----------------------------------------------------------------------------
## Estimating Area of Applicability (AOA)
## Area of Applicability

# We have seen that technically, the trained model can be applied to the entire area of interest 
# (and beyond...as long as the sentinel predictors are available which they are, even globally). 
# But we should assess if we SHOULD apply our model to the entire area. The model should only be 
# applied to locations that feature predictor properties that are comparable to those of the 
# training data. If dissimilarity to the training data is larger than the disimmilarity within 
# the training data, the model should not be applied to this location.
# 
# The calculation of the AOA is quite time consuming. To make a bit faster we use a parallelization.


cl <- makeCluster(4) # Using 4 clusters

registerDoParallel(cl)

AOA <- aoa(AsharVIstack, model, cl=cl)

plot(AOA)


## Plot and compare Prediction Error and AOA

# The result of the aoa function has two layers: the dissimilarity index (DI) and 
# the area of applicability (AOA). The DI can take values from 0 to Inf, where 0 means 
# that a location has predictor properties that are identical to properties observed in 
# the training data. With increasing values the dissimilarity increases. 
# The AOA has only two values: 0 and 1. 0 means that a location is outside the area of applicability, 
# 1 means that the model is inside the area of applicability. 
# Find more information on how the AOA is derived 
# in [Meyer\&Pebesma (2020)](http://arxiv.org/abs/2005.07939).



# predplot <- spplot(deratify(VIprediction), col.regions=cols, 
#                    main = list(label="Prediction (left), prediction only for the AOA (right) and RGB composite (bottom)",cex=0.8))

predplot <- spplot(deratify(VIprediction), col.regions=cols, 
                   main = list(label="Prediction (left) and prediction only for the AOA (right)",cex=0.8))

predplotaoa <- spplot(deratify(VIprediction), col.regions=cols)+
  spplot(AOA$AOA,col.regions=c("grey","transparent"))

# rgbplot_AsharVI <- spplot(AsharVIstack[[1]],  col.regions="transparent",
#                      sp.layout =rgb2spLayout(AsharVIstack[[1:3]], quantiles = c(0.02, 0.98), alpha = 1))

#latticeCombineGrid(list(predplot, predplotaoa, rgbplot_AsharVI),layout=c(2,2))
latticeCombineGrid(list(predplot, predplotaoa),layout=c(2,1)) # plot images in grid of 1row and 2columns

plot(predplotaoa)




library(sp)

grid.arrange(
  spplot(deratify(VIprediction),col.regions=cols, main="class prediction"),
  spplot(AOA$DI,col.regions=viridis(100), main="Dissimilarity Index"),
  spplot(deratify(VIprediction), col.regions=cols, main="Class prediction for AOA")+ 
    spplot(AOA$AOA,col.regions=c("grey","transparent")), ncol=3)


grid.arrange(spplot(AOA$DI,col.regions=viridis(100),main="Dissimilarity Index (DI)"),
             spplot(deratify(VIprediction), col.regions=cols, main="prediction for AOA \n(spatial CV error applies)")+
               spplot(deratify(AOA$AOA),col.regions=c("grey","transparent")), ncol=2)


grid.arrange(spplot(AOA$DI,col.regions=viridis(100),main="Dissimilarity Index (DI)"),
             spplot(deratify(VIprediction), col.regions=cols, main="Class prediction for AOA")+
               spplot(deratify(AOA$AOA),col.regions=c("grey","transparent")), ncol=2)

spplot(deratify(VIprediction), col.regions=cols, main="Class prediction for AOA")+
  spplot(deratify(AOA$AOA),col.regions=c("grey","transparent"))




## RESETTING THE tmap() properties
# reset the options to the default values
tmap_options_diff()

## tmap options successfully reset
tmap_options_reset()



#-------------------------------------------------------------------------------
#----------------------------------------
#Export ClusterMap as KML

# #Reproject the raster layer to WGS84 (<- "+init=epsg:28992") 				  
# # transform CRS to longitude/latitude
# # You need to reproject raster to lat/long to meet requirements for KML function
# Scenario1Cass_WGS <- projectRaster(Scenario1Cass, crs="+proj=longlat +datum=WGS84", method='ngb')
# 
# 
# # Export map using the KML() function.
# # Export raster data to a KML file and an accompanying PNG image file. 
# # Multi-layer objects can be used to create an animation. 
# # The function attempts to combine these into a single (and hence more convenient) 
# # KMZ file (a zip file containing the KML and PNG files).
# 


# Save Predicted Map in UTM Coordinates and as .tif file

writeRaster(VIprediction, file= paste0(Wdir, "/VIoutput/Scenario2allVIGLCM_MapUTM.tif"), datatype='FLT4S', overwrite=TRUE)

writeRaster(AOA$AOA, file= paste0(Wdir, "/VIoutput/Scenario2allVIGLCM_MapUTM_AOA.tif"), datatype='FLT4S', overwrite=TRUE)

writeRaster(AOA$DI, file= paste0(Wdir, "/VIoutput/Scenario2allVIGLCM_and_AOA_MapWGS.tif"), datatype='FLT4S', overwrite=TRUE)
AOA$parameters


VIGLCMvegMap <- raster(file.choose('Choose the classified VI-GLCM Scenario2 Map tiff file: '))



kmlColors <- c( "black", "turquoise", "khaki", "magenta", 
                "yellow", "purple", "cyan4", "sandybrown", "green", "sienna", "red")

# Setting the same colour codes as on map in QGIS

kmlColors2 <- c("#000201", "#72cea5", "#e69428", "#047c68", "#ec4f16", "#e8d22c", "#885bb5", 
                "#7d8b8f", "#3cec2c", "#91522d", "#d05e83")

KML(VIGLCMvegMap, file=paste0(Wdir, "/VIoutput/VIoutputScenario2/VIGLCMScenario2all_Prediction_AOA_MapWGS.kml"), col = kmlColors2, blur=1, zip='', overwrite=TRUE)









#--------------------------------------------------------
# Create Functions to Divide AOI into Distribution Grids
#--------------------------------------------------------

#------------------------------------------------------
# Function 1. To Divide Study Areas (or AOI) in Grids of a Specified Target Size

SpeciesGrids <- function(studyAOI, espgCode, gridSize){
  
  # Load the required libraries
  
  require(here) 
  require(rgdal) 
  require(tidyverse)
  require(raster)
  require(sf)
  require(sp)
  require(fasterize)
  require(rasterize)
  require(rgeos)
  require(dismo)
  
  require(dplyr)
  require(ggplot2)
  require(stringr) # use function to glue string
  
  studyAOI <- readOGR(file.choose("Select the shapefile for AlUla or Al-Madinah Province"))
  #studyAOI <- readOGR(file.choose("Select the shapefile for AOI"))
  espgCode <- as.numeric(readline(prompt = "Enter the ESPG Code for target CRS in UTM .e.g 32639 for Soqotra: "))
  gridSize <- as.numeric(readline(prompt = "Enter the target resolution or size for each grid e.g.2000: "))
  
  
  # Set the CRS for the AOI Shapefile
  
  #studyAOI <- spTransform(studyAOI, CRS(str_glue("+init=epsg:","{espgCode}")))
  studyAOI <- spTransform(studyAOI, CRS(str_glue("+init=epsg:{espgCode}")))
  
  
  # Create function for numbering of large grids
  #  to handle cases where row number exceeds 26
  
  excel_col <- function(n = NULL) {
    if(n < 1) { stop("Please provide a positive number")}
    if(n <= 26) { return(LETTERS[n]) } # return the corresponding UPPERCASE LETTER indexed by "n"
    out = c()
    while(n > 0) {
      rem <- n %% 26
      if(rem == 0) {
        out <-  append(out, 'Z')
        n <- (n/26) - 1
      } else {
        out <- append(out, LETTERS[rem])
        n <- n/26
      }
    }
    #paste0(rev(out), collapse = '')
    paste0(sort(out, decreasing = T), collapse = '')
  }
  
  
  
  # Make a grid across the whole AOI and the cell size that matches the target resolution
  
  ng2 <- sf::st_make_grid(studyAOI, cellsize = gridSize) 
  #ng2
  
  
  # Extract the bounding box from AOI shapefile
  reqGridBbox <- st_bbox(studyAOI)
  
  # Calculate number of rows/columns needed to hit your desired cellsize for a new grid (roundup values to nearest integer)
  nCols <- ceiling((reqGridBbox[3]-reqGridBbox[1])/gridSize) #xmax - xmin/gridSize
  nRows <- ceiling((reqGridBbox[4]-reqGridBbox[2])/gridSize) #ymax - ymin/gridSize
  
  
  # Create a sequence of ID for the length or row in the empty grid cell
  ng2 <- sf::st_sf(ng2, 'ID' = seq(length(ng2)), ng2)
  #ng2
  
  
  #Create new empty raster file
  # The raster file have to specify matching row/col and the ID of the empty grid
  # Specify the calucated rows and columns in both the grid and the raster conversion. 
  nr2 <- fasterize::raster(ng2, ncol = nCols, nrow = nRows) 
  nr2[] <- ng2$ID
  
  # Get X and Y from row and col numbers, convert Y to letters
  
  ng2$X <- raster::colFromCell(nr2, ng2$ID) # or seq(ncell(nr1))
  ng2$Y <- rev(sapply(raster::rowFromCell(nr2, ng2$ID), excel_col))
  
  # Combine the X and Y into new variable (LAB) as Label for each Grid
  ng2$LAB <- paste0(ng2$Y, ng2$X)
  
  
  
  st_crs(ng2) <- espgCode # set the CRS for the ng2
  #class(ng2)
  
  # convert sf multi-polygon data frame to large SpatialPolygonsDataFrame
  ng2_sp <- as(ng2, Class = "Spatial") 
  
  #str(ng2_sp)
  #plot(ng2_sp)
  
  
  # Set the projection CRS for the species occurrence grid to that of AOI shapefile
  projection(ng2_sp) <- projection(studyAOI)
  
  # intersect from raster package: 
  
  AOIintersectGrid <- intersect(ng2_sp, studyAOI) # clip the grid to boundary of AOI (Soqotra Mainland)
  AOIintersectGrid$Cellid <- 1:nrow(AOIintersectGrid) # Create Id values for all cells in intersect boundary
  
  ng2_sp$Cellid <- 1:nrow(ng2_sp) # Create Cell id for all square grids of the AOI extent
  
  # Export clipped grid as shapefile to specified directory
  dir.create("AOIShapefileGrid")
  
  writeOGR(obj <- ng2_sp, dsn = "AOIShapefileGrid", layer = "AOIUnclippedGrid", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
  writeOGR(obj <- AOIintersectGrid, dsn = "AOIShapefileGrid", layer = "AOIClippedGrid", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
  gridListAOI <- list("squareGrid" = ng2_sp, "clippedGrid" = AOIintersectGrid)
  
  return(gridListAOI)
  
}



## Apply Function 1 to a case: Example of Soqotra Archipelago - Main land Soqotra

SoqotraGrided <- SpeciesGrids() 

SoqotraGrided <- rgdal::readOGR(file.choose("Select the shapefile for AOI .i.e. AOIClippedGrid"))
plot(SoqotraGrided)
names(SoqotraGrided)







