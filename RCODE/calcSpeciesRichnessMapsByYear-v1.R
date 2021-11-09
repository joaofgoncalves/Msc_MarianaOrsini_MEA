

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

rm(list=ls())



source("./RCODE/_LibMainFunctions-v1.R")

# plotOutDir <- "./OUT/MODS/R02_Plots"
#plotOutDir <- "./OUT/MODS/R02_Plots_CropSE"
plotOutDir <- "./OUT/MODS/R02_Plots_CropRBSE"

# spGroups <- list.dirs("./OUT/MODS/R02", recursive = FALSE)
#spGroups <- list.dirs("./OUT/MODS/R02_CropSE", recursive = FALSE)
spGroups <- list.dirs("./OUT/MODS/R02_CropRBSE", recursive = FALSE)


#outDir <- "./OUT/MODS/R02_SpRichMapsByYear_cropSE"
outDir <- "./OUT/MODS/R02_SpRichMapsByYear_cropRBSE"

## ---------------------------------------------------------------------------------- ##

yrs_ <- seq(1985,2019,by=2)
yrs <- rep(seq(1985,2019,by=2),each=2)

selVarsDF <- readxl::read_excel("./DATA_/TABLES/LULC_Habitats_Species-simple-v2.xlsx")



for(spGroup in spGroups){
  k <- 0
  spNames <- list.files(spGroup, recursive = TRUE, full.names = TRUE, pattern=".tif$") 
  
  
  for(yr in yrs_){
    k <- k + 1
    spFilesBin <- spNames[grepl("TSSbin",spNames)]
    spFilesBinYr <- spFilesBin[grepl(paste("yr",yr,sep=""),spFilesBin)]
    
    spGroup_ <- basename(spGroup)
    
    
    rstBin <- rast(spFilesBinYr)
    rstSpRich <- app(rstBin, sum, na.rm=TRUE)
    
    fnOut <- paste(outDir,"/",spGroup_,"_spRichnessMap_",yr,"_v1.tif",sep="")
    
    if(k==1){
      allFiles <- fnOut
    }else{
      allFiles <- c(allFiles, fnOut)
    }
    
    writeRaster(rstSpRich, fnOut, overwrite = TRUE)
    
    cat("Finished map",spGroup_,"... year",yr,"\n\n")
  }
  
  rstAll <- rast(allFiles)
  rstSen <- app(rstAll, senTrendSlope, cores = 1)
 
  
  # fnOutSlope1 <- paste(outDir,"/_",spGroup_,"_spRichnessMap_SenSlope_cropSE_v1.tif",sep="") 
  # fnOutSlope2 <- paste(outDir,"/_",spGroup_,"_spRichnessMap_SenSlopePval_cropSE_v1.tif",sep="") 
  # fnOutSlope3 <- paste(outDir,"/_",spGroup_,"_spRichnessMap_SenSlopePvLTE005_cropSE_v1.tif",sep="") 
  
  fnOutSlope1 <- paste(outDir,"/_",spGroup_,"_spRichnessMap_SenSlope_cropRBSE_v1.tif",sep="")
  fnOutSlope2 <- paste(outDir,"/_",spGroup_,"_spRichnessMap_SenSlopePval_cropRBSE_v1.tif",sep="")
  fnOutSlope3 <- paste(outDir,"/_",spGroup_,"_spRichnessMap_SenSlopePvLTE005_cropRBSE_v1.tif",sep="")
  
  
  
  writeRaster(rstSen[[1]], fnOutSlope1, overwrite = TRUE)
  writeRaster(rstSen[[2]], fnOutSlope2, overwrite = TRUE)
  
  rstSen[[1]][rstSen[[2]] > 0.05] <- NA
  writeRaster(rstSen[[1]], fnOutSlope3, overwrite = TRUE)
  
}




