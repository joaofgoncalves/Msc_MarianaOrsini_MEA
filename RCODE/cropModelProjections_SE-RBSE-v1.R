

library(terra)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)


se <- vect("./DATA_/VECTOR/Area_estudo/espinhaco_total_Simplf_WGS84.shp")
rbse <- vect("./DATA_/VECTOR/2. Reserva_Biosfera_Serra_Espinhaco/Fase 1_2005/Limite_rbse_fase1/limite_dissolved_rbse1_buffer_30m.shp")


se_outDir <- "./OUT/MODS/R02_CropSE"
rbse_outDir <- "./OUT/MODS/R02_CropRBSE"



plotOutDir <- "./OUT/MODS/R02_Plots"


spGroups <- list.dirs("./OUT/MODS/R02", recursive = FALSE)

yrs <- seq(1985,2019,by=2)

selVarsDF <- readxl::read_excel("./DATA_/TABLES/LULC_Habitats_Species-simple-v2.xlsx")


for(spGroup in spGroups){
  
  spNames <- list.dirs(spGroup, recursive = FALSE) 
  
  spGroup_ <- basename(spGroup)
  
  outDir_SE <- paste(se_outDir,"/",spGroup_,sep="")
  outDir_RBSE <- paste(rbse_outDir,"/",spGroup_,sep="")
  
  if(!dir.exists(outDir_SE)){
    dir.create(outDir_SE)
  }
  if(!dir.exists(outDir_RBSE)){
    dir.create(outDir_RBSE)
  }
  
  
  for(spName in spNames){
    
    spName_ <- basename(spName)
    
    
    outDir_SE_spPath <- paste(se_outDir,"/",spGroup_,"/",spName_,sep="")
    outDir_RBSE_spPath <- paste(rbse_outDir,"/",spGroup_,"/",spName_,sep="")
    
    if(!dir.exists(outDir_SE_spPath)){
      dir.create(outDir_SE_spPath)
    }
    if(!dir.exists(outDir_RBSE_spPath)){
      dir.create(outDir_RBSE_spPath)
    }
    
    
    spNameFull <- selVarsDF$Species[selVarsDF$SpeciesComp %in% spName_]
    
    spFileList <- list.files(spName, pattern=".tif$", full.names = TRUE)
    
    binFiles <- spFileList[grepl("TSSbin",spFileList)] 
    hsbFiles <- spFileList[!grepl("TSSbin",spFileList)]
    
    rstBin <- rast(binFiles)
    rstHsb <- rast(hsbFiles)
    
    rstBin_cropSE <- mask(rstBin, se)
    rstBin_cropRBSE <- mask(rstBin, rbse)
    
    rstHsb_cropSE <- mask(rstHsb, se)
    rstHsb_cropRBSE <- mask(rstHsb, rbse)
    
    
    for(i in 1:nlyr(rstBin_cropSE)){
      
      fPath_BinSE <- paste(outDir_SE_spPath, basename(binFiles)[i],sep="/")
      fPath_BinRBSE <- paste(outDir_RBSE_spPath, basename(binFiles)[i],sep="/")
      
      fPath_HsbSE <- paste(outDir_SE_spPath, basename(hsbFiles)[i],sep="/")
      fPath_HsbRBSE <- paste(outDir_RBSE_spPath, basename(hsbFiles)[i],sep="/")
      
      writeRaster(rstBin_cropSE[[i]], fPath_BinSE, overwrite=TRUE)
      writeRaster(rstBin_cropRBSE[[i]], fPath_BinRBSE, overwrite=TRUE)
      
      writeRaster(rstHsb_cropSE[[i]], fPath_HsbSE, overwrite=TRUE)
      writeRaster(rstHsb_cropRBSE[[i]], fPath_HsbRBSE, overwrite=TRUE)
    }
    
    cat(":: Finished",spGroup_,"|",spNameFull,"::\n\n")
  }
  
  
}


