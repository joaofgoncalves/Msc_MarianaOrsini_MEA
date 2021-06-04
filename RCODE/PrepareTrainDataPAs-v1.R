
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(exactextractr)
library(stringi)
library(stringr)
library(rlang)
library(crayon)
library(randomForest)
library(biomod2)
library(biomod2plus)
library(trend)
library(diptest)

## -------------------------------------------------------------------------------------- ##

#setwd("C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA")

# Read ancillary functions
source("./RCODE/_LibMainFunctions-v1.R")

grid1k_wgs84 <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid_WGS84.shp")

dfPathList <- list.files("./DATA_/TABLES/_PRED_DATASETS/", ".rds$", full.names = TRUE)


## -------------------------------------------------------------------------------------- ##

spData <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Amphibians/endemic_amphibia_Espinhaco_v2-1.xlsx")


spDatasf <- prepSpData(spData, filterByMinObs = FALSE)

spDataGrid <- prepSpDataWithGridID(spDatasf, grid1k_wgs84,
                         filterByMinObs = FALSE, getCounts = TRUE)

spCountsAmph <- spDataGrid[["spCounts"]] %>% st_drop_geometry()
write_csv(spCountsAmph,"./OUT/Species1KmCounts_Amphibians-v1.csv")


createTrainData (spData           = spData,
                   spName         = NULL,
                   sfGrid         = grid1k_wgs84,
                   dfPathList     = dfPathList,
                   spNameCol      = "Species",
                   lonCol         = "Lon",
                   latCol         = "Lat",
                   yearCol        = "Year",
                   yrStart        = 1985,
                   yrEnd          = 2019,
                   filterByMinObs = TRUE,
                   nmin           = 20,
                   removeDups     = TRUE,
                   nPAsets        = 10,
                   nPAperSet      = "equal",
                   progressBar    = TRUE,
                   outDir         = "./DATA_/TABLES/_TRAIN_DATASETS/")

## -------------------------------------------------------------------------------------- ##

spData <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Birds/DH_Highland_Birds-v3.xlsx")


spDatasf <- prepSpData(spData, 
                       filterByMinObs = FALSE,
                       spNameCol      = "SpeciesName",
                       latCol         = "dLat",
                       lonCol         = "dLon",
                       yearCol        = "Year"
                       )

spDataGrid <- prepSpDataWithGridID(spDatasf, grid1k_wgs84,
                                   filterByMinObs = FALSE, getCounts = TRUE)

spCountsBirds <- spDataGrid[["spCounts"]] %>% st_drop_geometry()
write_csv(spCountsBirds,"./OUT/Species1KmCounts_HighlandBirds-v2.csv")


createTrainData (spData           = spData,
                 spName         = NULL,
                 sfGrid         = grid1k_wgs84,
                 dfPathList     = dfPathList,
                 spNameCol      = "SpeciesName",
                 lonCol         = "dLon",
                 latCol         = "dLat",
                 yearCol        = "Year",
                 yrStart        = 1985,
                 yrEnd          = 2019,
                 filterByMinObs = TRUE,
                 nmin           = 20,
                 removeDups     = TRUE,
                 nPAsets        = 10,
                 nPAperSet      = "equal",
                 progressBar    = TRUE,
                 outDir         = "./DATA_/TABLES/_TRAIN_DATASETS/")

