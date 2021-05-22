

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(fasterize)
library(exactextractr)
library(rasterDT)
library(readr)

#setwd("C:/MyFiles/Projects/Msc_MarianaOrsini_MEA")

setwd("C:/Users/JG/Google Drive/CIBIO_Desktop/R_projects/Msc_MarianaOrsini_MEA")

yrs <- 1985:2019

fl <- list.files("./DATA_/RASTER/LULC_1985_2019_int", 
                 pattern=".tif$", full.names = TRUE)

pb <- txtProgressBar(min=1, max=length(fl),style = 3)


for(i in 1:length(fl)){
  
  outGridFn <- paste("./DATA_/RASTER/RstGrids/rstGrid_LULC_",yrs[i],".tif",sep="")

  if(file.exists(outGridFn)){
    cat("\n\nSkipping file:\n-> ",outGridFn,"\n\n")
    next
  }else{
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
  }

  setTxtProgressBar(pb, i)
}




pb <- txtProgressBar(min=1, max=length(fl),style = 3)


for(i in 1:length(fl)){

    outFn <- paste("./DATA_/TABLES/VARS_PREP/LULC/CTab_LULC_",yrs[i],"_v1.csv" ,sep="")

    
  if(file.exists(outFn)){
    cat("\n\nSkipping file:\n-> ",outFn,"\n\n")
    next
  }else{
    outGridFn <- paste("./DATA_/RASTER/RstGrids/rstGrid_LULC_",yrs[i],".tif",sep="")
    rstGrid <- raster(outGridFn)
    rstLULC <- raster(fl[i])

    rstCtab <- rasterDT::crosstabDT(rstGrid, rstLULC)
    cnames <- intersect(as.character(1:17), colnames(rstCtab))
    rstCtab <- rstCtab[,cnames]
    
    colnames(rstCtab) <- paste("cl",colnames(rstCtab),sep="_")
    rstCtab <- data.frame(ID = 1:nr, rstCtab)
    readr::write_csv(rstCtab, outFn)
  }

  setTxtProgressBar(pb, i)
  
}

