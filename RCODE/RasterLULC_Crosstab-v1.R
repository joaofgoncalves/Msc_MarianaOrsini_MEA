

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(fasterize)
library(exactextractr)
library(rasterDT)
library(readr)

setwd("C:/MyFiles/Projects/Msc_MarianaOrsini_MEA")

yrs <- 1985:2019

fl <- list.files("./DATA_/RASTER/LULC_1985_2019_int", 
                 pattern=".tif$", full.names = TRUE)

pb <- txtProgressBar(min=1, max=length(fl),style = 3)


for(i in 1:length(fl)){
  
 
  if(exists(outGridFn)){
    cat("\n\nSkipping file:\n-> ",outGridFn,"\n\n")
    next
  }
  
  outGridFn <- paste("./DATA_/RASTER/RstGrids/rstGrid_LULC_",yrs[i],".tif",sep="")

  
  grd <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid_WGS84.shp")
  nr <- nrow(grd)
  
  rstLULC <- raster(fl[i])
  rstGrid <- fasterize(grd, rstLULC, field="ID")
  rm(grd)
  
  writeRaster(rstGrid, 
              filename  = outGridFn, 
              overwrite = TRUE, 
              datatype  = "INT4S")
  
  rm(rstGrid)
  setTxtProgressBar(pb, i)
}




pb <- txtProgressBar(min=1, max=length(fl),style = 3)


for(i in 1:length(fl)){

  outGridFn <- paste("./DATA_/RASTER/RstGrids/rstGrid_LULC_",yrs[i],".tif",sep="")
  
  rm(rstGrid)
  rstGrid <- raster(outGridFn)
  
  rstCtab <- rasterDT::crosstabDT(rstGrid, rstLULC)
  cnames <- intersect(as.character(1:17), colnames(rstCtab))
  rstCtab <- rstCtab[,cnames]
  
  colnames(rstCtab) <- paste("cl",colnames(rstCtab),sep="_")
  rstCtab <- data.frame(ID = 1:nr, rstCtab)
  
  outFn <- paste("./DATA_/TABLES/VARS_PREP/LULC/CTab_LULC_",yrs[i],"_v1.csv" ,sep="")
  
  readr::write_csv(rstCtab, outFn)
  
  setTxtProgressBar(pb, i)
  
}

