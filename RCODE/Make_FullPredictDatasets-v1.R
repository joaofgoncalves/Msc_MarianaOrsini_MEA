

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(exactextractr)
library(stringi)
library(stringr)
library(rlang)

staticDataPaths <- matrix(
  c("clim","./DATA_/TABLES/VARS/BioClimDF-v1.rds",
    "topo","./DATA_/TABLES/VARS/TopoDF-v1.rds",
    "hydr","./DATA_/TABLES/VARS/HydroDF-v1.rds"),
  byrow = TRUE, ncol=2)

dynamicDataPaths <- matrix(
  c("EVI2","./DATA_/TABLES/VARS/EVI2DF_ByYear1985_2019-v1.rds",
   "NDWI","./DATA_/TABLES/VARS/NDWIDF_ByYear1985_2019-v1.rds",
   "LULC","./DATA_/TABLES/VARS/LULCDF_ByYear1985_2019-v1.rds"),
  byrow = TRUE, ncol=2)

outDir <- "./DATA_/TABLES/_PRED_DATASETS/"


clim <- read_rds(staticDataPaths[1,2])
topo <- read_rds(staticDataPaths[2,2])
hydr <- read_rds(staticDataPaths[3,2])

staticDF <- clim %>% 
  left_join(topo, by="ID") %>% 
  left_join(hydr, by="ID")

rm(list = c("clim","topo","hydr"))

EVI2 <- read_rds(dynamicDataPaths[1,2])
NDWI <- read_rds(dynamicDataPaths[2,2])
LULC <- read_rds(dynamicDataPaths[3,2])


yrs <- 1985:2019

pb <- txtProgressBar(1,length(yrs),style=3)

for(i in 1:length(yrs)){
  
  suffName <- paste("_",yrs[i],sep="")
  
  dynamicTMP_DF <- staticDF %>% 
    left_join((EVI2 %>% select(ID, ends_with(suffName))), by="ID") %>% 
    left_join((NDWI %>% select(ID, ends_with(suffName))), by="ID") %>% 
    left_join((LULC %>% select(ID, ends_with(suffName))), by="ID")
  
  newColNames <- gsub(suffName,"",colnames(dynamicTMP_DF))
  colnames(dynamicTMP_DF) <- newColNames
  outPath <- paste(outDir,"PredDF_SEspinhacoBBoxGrid1K_",yrs[i],".rds",sep="")
  
  write_rds(dynamicTMP_DF,outPath)
  setTxtProgressBar(pb,i)
}
