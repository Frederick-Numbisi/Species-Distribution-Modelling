
# SOQOTRA SPECIES DISTRIBUTION DATA MINING AND GRID CLASSIFICATION

#--------------------------------------------------------
# Create Functions to Classify Species Distribution Grids
#--------------------------------------------------------

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

setwd('/Users/fnumbisi/Documents/RBGE_MajorFloras_CMEP/SOQOTRA_Project/')


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
  
  studyAOI <- readOGR(file.choose("Select the shapefile for AOI"))
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

SoqotraGrided <- readOGR(file.choose("Select the shapefile for AOI .i.e. AOIClippedGrid"))
plot(SoqotraGrided)
names(SoqotraGrided)





#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# FUNCTION 2. CLASSIFYING SPECIE DISTRIBUTION GRIDS
# 1. Extract and create respective lists of shapefiles of species' Inferred and Observed Distributions
# 2. Estimating Area for each grid of the Inferred distribution
# 3. Use Observed distribution to classify Grids of Inferred Distribution into 1 of 3 Classes:
#    Actual Presence, Inferred Presence, and Absence.

#-------------------------------------------------------------------------------------



SoqWdir <- "/Users/fnumbisi/Documents/RBGE_MajorFloras_CMEP/SOQOTRA_Project/"


# Extract and create list of observed and inferred distribution
# Extract all the observed distribution files in a directory with .shp (shapefile) extension


# --------------------------------------------------
# Extract files with particular (desired) extensions
# Note: use $ at the end implies this is end of string to ensure that you match only files with desired extension
# For example, "\\.shp$" will exclude files ending with .shp.xml, which are created in QGIS and ArcGIS table of polygons
# Similarly, .dbf "dbf$" will only select geodatabase file and exclude the additional .dbf.xml files created in ArcGIS
# Also, you may use \\ to excape the special character . (dot). This is, however, optional.

InferrDistList <- sort(list.files(path = paste0(SoqWdir, "Soqotra_KBAs20220129/InferredDistPolygons"),full.names=TRUE, pattern='\\.shp$'))
ObserDistList <- sort(list.files(path = paste0(SoqWdir, "Soqotra_KBAs20220129/ObservedDistributionFiles"),full.names=TRUE, pattern='\\.shp$'))


# ------------------------------------------------
## Create Function for Pixel (Grid) Classification

SpeciesGridClassification <- function(LabeledGrid, espgCode, InferrDistList, ObserDistList, threshold){
  
  # Import the Labelled Clipped Grid for the AOI (Soqotra) Clipped with Cellid column
  LabeledGrid <- readOGR(file.choose("Select the Clipped AOI Grid Shapefile with 'Cellid' attribute i.e. AOIClippedGrid.shp"))
  
  # Get the ESPG Code for UTM CRS for the AOI
  crsCode <- as.integer(readline(prompt = "Enter the ESPG Code for target CRS in UTM .e.g 32639 for Soqotra: "))
  espgCode <- as.numeric(crsCode)
  
  LabeledGrid <- spTransform(LabeledGrid, CRS(str_glue("+init=epsg:{espgCode}")))
  
  threshold <- as.integer(readline(prompt = "Enter threshold percentage for area absence polygon.e.g 25 for 25% threshold: "))
  
  ## Create Directory to Save Results
  
  if (!dir.exists("SOQOTRAgridClass")) dir.create("SOQOTRAgridClass", recursive = TRUE)
  
  #dir.create("SOQOTRAgridClass")

  
  ## Create Duplicate of AOI (Soqotra) Grid
  
  LabGrid <- LabeledGrid
  
  FinalLabGrid <- c()
  
 
  outPut <- list()
  
  SpeciesList <- list()
 
  # --------------------------------
  # Use a sing placeholder variable, i, to iterate over both lists (Inferred and Observed Distribution)
  # Iterate through all the files of Inferred and Observed Species Distribution
  
  for (i in 1:length(InferrDistList)){
      
    #results <- vector(mode = "integer", length = length(ObserDistList))
      #data <- read.csv(InferrDistList[i]) # read the ith file from the list by indexing
      #countRows <- nrow(data) # count the number of rows in the file
      #results[i] <- countRows # store the counted rows in the ith position of the results vector
      
    
      # --------------------------------
      # Import the Inferred and Observed Species distribution
      
      #infDist <- readOGR(file.choose("Select the shapefile for Inferred Species Distribution: "), stringsAsFactors = F)
      
      infDist <- readOGR(InferrDistList[i], stringsAsFactors = F)
      #infDist <- readOGR(InferrDistList[i])
      
      
      SpeciesA_inf<- spTransform(infDist, CRS(str_glue("+init=epsg:{espgCode}")))
      
      # Get just the file name without the extension
      # remove the .shp extension using the str_remove()
      # particular shapefile pattern
      NameSpecies <- str_remove(basename(InferrDistList[i]), pattern = ".shp") 
      
      
      
      #plantSpecies <- basename(file.choose("choose the observed distribution points shapefile: "))
      
      obsDist <- readOGR(ObserDistList[i], stringsAsFactors = F)
      #obsDist <- readOGR(ObserDistList[j])
      
      SpeciesAObs<- spTransform(obsDist, CRS(str_glue("+init=epsg:{espgCode}")))
      
      #speciesName <- str_remove(plantSpecies, pattern = '.shp') # extract name of species
      
      
      # Create the %notin% operator since it is not built operator in R
      `%notin%` <- Negate(`%in%`) 
      
      
      # A) Create an intersect of SoqotraClippedGrid with the inferred and observed Species Distribution
      
      # 1) Intersect with clipping to polygon boundary for inferred distribution
      
      Clippedgrid_SpeciesAinf <- intersect(SpeciesA_inf, LabeledGrid)
      plot(Clippedgrid_SpeciesAinf)
      
      # 2) Intersect without clipping for inferred distribution: retain square grid to intersect cells
      
      SpeciesInfGrid <- raster::intersect( LabeledGrid, SpeciesA_inf)
      plot(SpeciesInfGrid)
      
      grid_SpeciesAinf <- LabeledGrid[LabeledGrid$Cellid %in% SpeciesInfGrid$Cellid, ]
      #plot(grid_SpeciesAinf)
      #grid_SpeciesAinf$Cellid
      
      # 3) Intersect without polygon clip for observed distribution
      #SpeciesObsGrid <- raster::intersect(SpeciesAObs, LabeledGrid) # Observation points intersect grid
      SpeciesObsGrid <- intersect(LabeledGrid, SpeciesAObs) # Observation points intersect grid
      #plot(SpeciesObsGrid)
      
      grid_SpeciesAobs <- LabeledGrid[LabeledGrid$Cellid %in% SpeciesObsGrid$Cellid, ]
      #plot(grid_SpeciesAobs)
      #grid_SpeciesAobs$Cellid
      
      # B) Estimate the total area for each polygon of the clipped inferred distribution grid
      
      # Estimate the area and save in new variable "Area_m2"
      # Extract areas from polygon objects then attach as attribute
      
      Clippedgrid_SpeciesAinf$Area_m2 <- area(Clippedgrid_SpeciesAinf)
      
      # Filter the Grids that fall in the target pixel classes "Absence" 
      # Use the estimated area to filter grids in which the inferred polygon covers < 25% of total grid size (4 Squared Km)
      
      Threshold = (threshold/100)*4e+06
      ClippedGridClass3 <- Clippedgrid_SpeciesAinf[Clippedgrid_SpeciesAinf$Area_m2 < Threshold & Clippedgrid_SpeciesAinf$Cellid %notin%grid_SpeciesAobs$Cellid, ]
      #ClippedGridClass3 <- Clippedgrid_SpeciesAinf[Clippedgrid_SpeciesAinf$Area_m2 < 1e+06, ]
      #plot(ClippedGridClass3)
      #ClippedGridClass3$Cellid
      
      
      
      # Use created operator (%notin%) to filter grids for the pixel class "Inferred Presence"
      
      #ClippedGridClass2 <- grid_SpeciesAinf[grid_SpeciesAinf$Cellid %notin%grid_SpeciesAobs$Cellid & !(grid_SpeciesAinf$Cellid %in% ClippedGridClass3$Cellid),  ]
      ClippedGridClass2 <- grid_SpeciesAinf[grid_SpeciesAinf$Cellid %notin%grid_SpeciesAobs$Cellid & grid_SpeciesAinf$Cellid %notin% ClippedGridClass3$Cellid,  ]
      
      #plot(ClippedGridClass2)
      #ClippedGridClass2$Cellid
      
      
      #------------------------------------------------------
      # Extract the Cell Ids for the respective pixel classes
      gridClass1id <- grid_SpeciesAobs$Cellid
      
      gridClass3id <- ClippedGridClass3$Cellid
      
      gridClass2id <- ClippedGridClass2$Cellid
      
      #gridClass1id%in%gridClass2id
      #gridClass1id%in%gridClass3id
      #gridClass2id%in%gridClass3id
      
      # ---------------------------------------------------------
      # Create a New Dataframe by extracting the Classid column
      # rm(Squaregrid_SpeciesAclassID) # remove df from workspace
      grid_SpeciesAclassID <- data.frame(grid_SpeciesAinf$Cellid)
      
      # use the dplyr package to rename the column
      #  as there is only one column in the New Dataframe, simply rename columns to desired name
      # alternatively use: names(df)[names(df) == "old.var.name"] <- "new.var.name"
      
      colnames(grid_SpeciesAclassID) <- "Cellid"
      
      
      # Create a new variable (PixelClass) in the New Dataframe using the %in% operator
      grid_SpeciesAclassID <- within(grid_SpeciesAclassID, {
        pixelClass = "Absence"
        pixelClass[Cellid %in% gridClass1id] = "Actual Presence"
        pixelClass[Cellid %in% gridClass2id] = "Inferred Presence"
        pixelClass[Cellid %in% gridClass3id] = "Absence"
      })
      
      # Now merge the squaredgrid for speciesA with the New Dataframe containing "pixelClass"
      #  merge by the attribute "Cellid" while preserving the class as SpatialPolygonDataFrame
      
      grid_SpeciesAinf <- merge(grid_SpeciesAinf, grid_SpeciesAclassID, by = "Cellid")
      
      # Extend the New Dataframe by adding attribute Geometry to the dataframe 
      # Add geometry from the Inferred Distribution Squared grid and match the polygons by Cellid
      # spdfSpeciesAPixels <- addAttrToGeom(grid_SpeciesAinf, grid_SpeciesAclassID, match.ID = "Cellid") 
      
      
      # Export spatialPolygonDataframe with pixel classes of species distribution as shapefile
      # Create directory to save results 
      
      # writeOGR(obj = grid_SpeciesAinf, dsn = "SoquotraGrid", layer = paste(NameSpecies,'PixelsClasses', sep = "_"), 
      #         driver = "ESRI Shapefile", overwrite_layer = TRUE )
      
      writeOGR(obj = grid_SpeciesAinf, dsn = "SoquotraGrid", layer = paste0(c(NameSpecies,'PixelsClasses'), collapse = "_"), 
               driver = "ESRI Shapefile", overwrite_layer = TRUE )
      
      
      gridClasses <- grid_SpeciesAinf
      # gridClasses <- list(grid_SpeciesAinf, spdfSpeciesAPixels)
      
      outPut <- append(outPut, gridClasses)
      
      SpeciesList <- append(SpeciesList, NameSpecies)
      
      # -----------------------------------------------------------------------
      
      
      ## Create new attribute column for Species and assign NA values
      LabGrid$Species <- NA
      
      
      ## Replace Species (NA) values with pixelClass for matching "Cellid"
      
      # Example df1$status[df1$ID %in% df2$ID] <- df2$status[df2$ID %in% df1$ID]
      
      LabGrid$Species[LabGrid$Cellid %in% grid_SpeciesAinf$Cellid] <- grid_SpeciesAinf$pixelClass[grid_SpeciesAinf$Cellid %in% LabGrid$Cellid] 
      
      
      ## To Coerce the sp object to sf object.
      # x <- as(x, "sf") 
      
      LabGrid2 <- as(LabGrid, "sf")
      
      # LabeledGrid_update1 <- as(LabeledGrid_update1, "sf")
      
      ## Rename the Species attribute to the corresponding Name of Species (using the variable "NameSpecies")
      # With the newest dplyr (>0.7.0) you would use
      # The syntax is "new_name=old_name" and you need to use := with !! 
      #  to put a variable on the left side of a parameter name.
      
      LabGrid2 <- rename(LabGrid2, !!NameSpecies:=Species)
      
      
      ## Coerce into spatial data frame as it is currently a SF
      # newPoly <- as(newData, "Spatial")
      
      LabGrid3 <- as(LabGrid2, "Spatial")
      
      
      LabGrid <- LabGrid3
      
      FinalLabGrid <- LabGrid2
      
      
       
    }
  
  
  # Get the column index (number) from the column name (in a dataframe)
  # which(colnames(df)=="b" )
  
  
  
  #writeOGR(obj = FinalLabGrid, dsn = "SOQOTRAgridClass", layer = paste0(c('sfAllSpeciesGridClasses'), collapse = "_"), 
  #         driver = "ESRI Shapefile", overwrite_layer = TRUE )
  
  writeOGR(obj = LabGrid, dsn = "SOQOTRAgridClass", layer = paste0(c('spAllSpeciesGridClasses'), collapse = "_"), 
           driver = "ESRI Shapefile", overwrite_layer = TRUE )
  
  write.table(FinalLabGrid, file= paste0(SoqWdir, "/SOQOTRAgridClass/sfAllSpeciesGridClasses.csv"), 
              append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
  
  return(list(outPut, FinalLabGrid, LabGrid))
}


## -----------------------------------------------------------------------------
## Run Classification Function for the Sample Species Distribution Lists

TrialSpeciesDist <- SpeciesGridClassification(InferrDistList=InferrDistList, 
                                              ObserDistList=ObserDistList)


## View the Combined Dataframe of Labled Grids for all species analysis
View(TrialSpeciesDist[2][[1]])

plot(TrialSpeciesDist[3][[1]])




