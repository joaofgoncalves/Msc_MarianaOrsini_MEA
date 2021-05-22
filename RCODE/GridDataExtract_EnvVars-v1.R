
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(exactextractr)
library(stringi)
library(stringr)
library(rlang)


## ---------------------------------------------------------------------- ##
## DATA
## ---------------------------------------------------------------------- ##


grid1k <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid.shp")


grid1k_wgs84 <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid_WGS84.shp")


TWI <- raster("./DATA_/RASTER/TWI_srtmv4.tif")

mdt <- raster("./DATA_/RASTER/mdt_bbox_utm.tif")
DNI <- raster("./DATA_/RASTER/Brazil_GlobalSolarAtlas_v2/DNI.tif")
DFI <- raster("./DATA_/RASTER/Brazil_GlobalSolarAtlas_v2/DIF.tif")



## ---------------------------------------------------------------------- ##
## PROCESS DATA
## ---------------------------------------------------------------------- ##


TWI_avg <- exact_extract(x   = TWI,
                         y   = grid1k,
                         fun = "mean")
TWI_std <- exact_extract(x   = TWI,
                         y   = grid1k,
                         fun = "stdev")

HydroDF <- 
  bind_cols(grid1k %>% st_drop_geometry %>% select(ID),
            TWI_avg = TWI_avg,
            TWI_std = TWI_std )

write_csv(HydroDF,"./DATA_/TABLES/VARS/HydroDF-v1.csv")
write_rds(HydroDF,"./DATA_/TABLES/VARS/HydroDF-v1.rds")

## ---------------------------------------------------------------------- ##

mdt_avg <- exact_extract(x   = mdt,
                         y   = grid1k,
                         fun = "mean")

mdt_std <- exact_extract(x   = mdt,
                         y   = grid1k,
                         fun = "stdev")

DNI_avg <- exact_extract(x   = DNI,
                         y   = grid1k_wgs84,
                         fun = "mean")

DFI_avg <- exact_extract(x   = DFI,
                         y   = grid1k_wgs84,
                         fun = "mean")

TopoDF <- 
bind_cols(grid1k %>% st_drop_geometry %>% select(ID),
          Elev_Avg = mdt_avg,
          Elev_Std = mdt_std,
          DNI_avg  = DNI_avg,
          DFI_avg  = DFI_avg)

write_csv(TopoDF,"./DATA_/TABLES/VARS/TopoDF-v1.csv")
write_rds(TopoDF,"./DATA_/TABLES/VARS/TopoDF-v1.rds")


## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##

#library(parallel)
#library(doParallel)
#registerDoParallel(4)


fl <- list.files("./DATA_/RASTER/CLIM_CHELSA/", pattern=".tif$", full.names = TRUE)

#foreach(i = 1:19, .combine = "cbind"){
for(i in 1:19){
    
  bioRst <- raster(fl[i])
  
  tmpDF <- exact_extract(x   = bioRst,
                         y   = grid1k_wgs84,
                         fun = "mean")
  
  vn <- paste("BIO",str_pad(i,2,"left","0"),sep="_")
  
  if(i==1){
    BioClimDF <- 
      bind_cols(grid1k %>% st_drop_geometry %>% select(ID),
      !!vn := tmpDF)
  }else{
    BioClimDF <- 
      bind_cols(BioClimDF, !!vn := tmpDF)
  }
}


write_csv(BioClimDF,"./DATA_/TABLES/VARS/BioClimDF-v1.csv")
write_rds(BioClimDF,"./DATA_/TABLES/VARS/BioClimDF-v1.rds")


