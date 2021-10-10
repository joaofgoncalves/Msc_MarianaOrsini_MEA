

library(readxl)
library(dplyr)
library(raster)
library(terra)
library(sf)

tr <- read_excel("./DATA_/TABLES/transitions_map.xlsx",col_names = FALSE) %>% as.matrix()

from <- tr[-1,1]
to <- tr[1,-1]
from==to

tr <- tr[-1,-1]

fromMap <- raster("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_1985_int.tif")
toMap <- raster("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_2019_int.tif")

outMap <- fromMap
#outMap[!is.na(outMap)] <- 0

outTransMatrix <- matrix(NA,ncol=3,nrow=length(from)^2)

pb <- txtProgressBar(1,length(from)^2,style=3)

k <- 0

for(i in 1:length(from)){
  
  fromClass <- as.integer(from[i])
  
  for(j in 1:length(to)){
    
    k <- k + 1
    
    toClass <- as.integer(to[j])
    outValue <- as.integer(tr[i,j])
      
    outTransMatrix[k,1] <- fromClass
    outTransMatrix[k,2] <- toClass
    outTransMatrix[k,3] <- outValue
    
    #fromMap == fromClass
    # print(fromClass)
    # print(toClass)
    # print(outValue)
    # print("-------------------")
    
    # mask1 <- (fromMap == fromClass)
    # mask2 <- (toMap == toClass)
    # 
    # outMap[mask1 & mask2] <- outValue
    setTxtProgressBar(pb, k)
  }
  
}

TRANS_MATRIX <<- outTransMatrix

TRANS_VEC <<- TRANS_MATRIX[,3]
names(TRANS_VEC) <- paste(TRANS_MATRIX[,1],TRANS_MATRIX[,2],sep="_")


# transFunc <- function(x,...){
#   #print(x)
#   #print(y)
#   ind <- try((TRANS_MATRIX[,1] == x) & (TRANS_MATRIX[,2] == y))
#   if(inherits(ind,"try-error")){
#     return(NA)
#   }else{
#     return(as.numeric(TRANS_MATRIX[ind,3]))
#   }
# } 
# 
# 
# transFunc <- function(x,...){
#   TRANS_VEC[paste(x,sep="_")]
# }
# 
# 
# 
# transRst <- overlay(fromMap, toMap, fun = transFunc)




transFunc <- function(x, y, ...){
  return(TRANS_VEC[paste(x,y,sep="_")])
}

fromMap <- rast("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_1985_int.tif")
toMap <- rast("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_2019_int.tif")

twoMaps <- c(fromMap, toMap)

transMap <- lapp(twoMaps, transFunc)

fromMap <- rast("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_2000_int.tif")

twoMaps <- c(fromMap, toMap)

transMap2000 <- lapp(twoMaps, transFunc)


terra::writeRaster(transMap,"./DATA_/RASTER/TransitionsMap_1985_2019_v1.tif", datatype="INT2S")


terra::writeRaster(transMap2000,"./DATA_/RASTER/TransitionsMap_2000_2019_v1.tif", datatype="INT2S")




