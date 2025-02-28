
library(basemaps) # to load basemaps
library(leaflet)
library(dismo) # also loads the sp and raster packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(raster)
library(sp)

theme_set(theme_minimal())

library(rgdal)
library(sf)
library(maptools) # to add north arrow
library(leaflet.extras)


setwd("/home/numbisi/Documents/pathToWorkingDirectory")

filepath1 = '/home/numbisi/Documents/filePath1'
filepath2 = '/home/numbisi/Documents/filePath2'


# Create function to exclude row with duplicate coordinates
#DupCoordinates = function(SpeciesData){
#  for(i in 3:ncol(SpeciesData)){
#    if (!is.na(SpeciesData[, i])){
#      SpeciesData[, i] = 1
#    } else {
#      SpeciesData[, i] = 0
#    }
#  }
#}


## Create a function to clean and recode the species occurrence data
### Species occurrence converted to binary: 0 (Absence) and 1 (Presence)
Occurrence2 = function(SpeciesData){
  for(i in 3:ncol(SpeciesData)){
    SpeciesData[, i][SpeciesData[, i] == "" | SpeciesData[, i] == " "] <- 0
    SpeciesData[, i][SpeciesData[, i] != 0] <- 1
  }
  return(SpeciesData)
}

# Create function to convert CRS from Degree Minute Seconds (DMS) to Decimal Degrees 
Degree2dec <- function(angle) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, split=' '))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    if (length(y) ==1) {
      y[1]  # for angle coordinates in Decimal Degrees
    } else if (length(y) ==2) {
      y[1] + y[2]/60 # for angle coordinates in Degree Decimal Minutes
    } else {
      y[1] + y[2]/60 + y[3]/3600 # computation for angle coordinates in DMS
    }
    
  })
  return(x)
}



psbrsShape = readOGR("./PMBSR/PMBSRDA.shp") # import shapefile of project extent


SpeciesMap = function(data1, startcolumn){
  
  # for-loop over columns beginning from first species
  psbrsShape = readOGR(file.choose("Select the shapefile for AOI")) # Prompt user to select and Load the shapefile for study area from relative directory
  
  data1 = data.frame(read.csv(file.choose('Chose the species occurrence csv file: ')))
  
  startcolumn = as.integer(readline(prompt="Enter the index number of first species colum: ")) # prompt user to input in number
  
  
  ## Create a function to clean and recode the species occurrence data
  
  Occurrence2 = function(SpeciesData){
    for(i in startcolumn:ncol(SpeciesData)){
      SpeciesData[, i][SpeciesData[, i] == "" | SpeciesData[, i] == " "] <- 0 # Allocate value of zero to empty cells
      SpeciesData[, i][SpeciesData[, i] != 0] <- 1 # Allocate value of 1 to any cell with species record ()
    }
    return(SpeciesData)
  }
  
  data1 = Occurrence2(data1) # Run function of species occurrrence data
  
  # Iterate over the species columns to subset the dataset be species
  for(i in startcolumn:ncol(data1)) {   
    dataSp = data1[data1[[i]] == 1,] # filter and subset rows for values of i column = 1
    Name = names(data1[i]) # return the name of column i
    dataSp = dataSp[, c('Latitude', 'Longitude', Name)] # subset data to contain 3 columns: Lat, Long and column i
    
    # Create directory for each species if necessary (would raise warning if directory already exists)
    #
    if (!dir.exists(Name)) dir.create(Name, recursive = TRUE)
    
    # Finally save each species data as csv within the species directory with files named according to the column names
    write.csv(dataSp, file.path(Name, paste0(Name, ".csv")), row.names = F) # Exclude rownames in csv files
    
    # Create map in leaflet
    # Call the required libraries
    require(rgdal)
    require(leaflet)
    require(mapview)
    require(dplyr)
    require(leaflet.extras)
    
    # Convert Latitude and Longitude to Numeric variables
    #transform(dataSp, Longitude = as.numeric(Longitude),
    #         Latitude = as.numeric(Latitude))
    
    dataSp[, 1:3] <- sapply(dataSp[, 1:3], as.numeric)
    
    m = leaflet(dataSp) %>%
      addProviderTiles(providers$Esri.WorldImagery ) %>%
      addPolygons(data = psbrsShape, color = "#000000", weight = 2, opacity = 1, fillOpacity = 0)%>%
      addCircleMarkers(data = dataSp, lat = ~Latitude, lng = ~Longitude, color = "#FF0000 ", fill = T, fillColor ="#FF0000 ", fillOpacity = 1,  radius = ~2)%>%
      addScaleBar(position = "bottomright", options = scaleBarOptions(maxWidth = 200, metric = T))
    
    # Save the Species Distribution Map in both PNG and HTML (Interactive) formats
    mapshot(m, file = paste0(getwd(), "_", Name,".png"), url = paste0(getwd(), "/", Name,".html"),
            remove_controls = c("homeButton", "layersControl", "zoomControl"))
  }
  
}


# Run the Species Distribution Map Function for the species occurence table, specifying the start column

SpeciesMap(PMBSRSpeciesUpdate, startcolumn)





