
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
library(terra)


## -------------------------------------------------------------------------------------- ##

#setwd("C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA")

# Read ancillary functions
source("./RCODE/_LibMainFunctions-v1.R")

#grid1k_wgs84 <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid_WGS84.shp")

grid1k_wgs84 <- vect("./DATA_/VECTOR/Area_estudo/espinhaco_grid_WGS84.shp")

dfPathList <- list.files("./DATA_/TABLES/_PRED_DATASETS/", ".rds$", full.names = TRUE)


## -------------------------------------------------------------------------------------- ##

spData <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Amphibians/endemic_amphibia_Espinhaco_v2-1.xlsx")


spDatasf <- prepSpData(spData, filterByMinObs = FALSE)

spDataGrid <- prepSpDataWithGridID(spDatasf, grid1k_wgs84,
                         filterByMinObs = FALSE, getCounts = TRUE)

spCountsAmph <- spDataGrid[["spCounts"]]

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
                   nmin           = 15,
                   removeDups     = TRUE,
                   nPAsets        = 10,
                   nPAperSet      = "equal",
                   progressBar    = TRUE,
                   outDir         = "./DATA_/TABLES/_TRAIN_DATASETS/v3/AMPHIBIANS")

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

spCounts <- spDataGrid[["spCounts"]]

write_csv(spCounts,"./OUT/Species1KmCounts_HighlandBirds-v2.csv")


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
                 nmin           = 15,
                 removeDups     = TRUE,
                 nPAsets        = 10,
                 nPAperSet      = "equal",
                 progressBar    = TRUE,
                 outDir         = "./DATA_/TABLES/_TRAIN_DATASETS/v3/BIRDS")


## -------------------------------------------------------------------------------------- ##

spData <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Reptiles/Repteis-Espinhaco_Database_2021-vii-15-Henrique_MO_JG_v3.xlsx")


spDatasf <- prepSpData(spData,
                       filterByMinObs = FALSE,
                       spNameCol      = "Species",
                       latCol         = "Latitude",
                       lonCol         = "Longitude",
                       yearCol        = "Year"
)

spDataGrid <- prepSpDataWithGridID(spDatasf, grid1k_wgs84,
                                   filterByMinObs = FALSE, getCounts = TRUE)

spCounts <- spDataGrid[["spCounts"]]

write_csv(spCounts,"./OUT/Species1KmCounts_Reptiles-v1.csv")


createTrainData (spData         = spData,
                 spName         = NULL,
                 sfGrid         = grid1k_wgs84,
                 dfPathList     = dfPathList,
                 spNameCol      = "Species",
                 lonCol         = "Longitude",
                 latCol         = "Latitude",
                 yearCol        = "Year",
                 yrStart        = 1985,
                 yrEnd          = 2019,
                 filterByMinObs = TRUE,
                 nmin           = 15,
                 removeDups     = TRUE,
                 nPAsets        = 10,
                 nPAperSet      = "equal",
                 progressBar    = TRUE,
                 outDir         = "./DATA_/TABLES/_TRAIN_DATASETS/v2/REPTILES")



