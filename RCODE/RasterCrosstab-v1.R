

library(raster)
library(sf)

r1 <- raster("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_1985_int.tif")

r2 <- raster("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_2019_int.tif")

esp <- read_sf("C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA/DATA_/VECTOR/Area_estudo/espinhaco_total_Simplf_WGS84.shp")


changeMatrix <- crosstab(r1, r2)


r1_Esp <- mask(r1, esp)

r2_Esp <- mask(r2, esp)


changeMatrix <- crosstab(r1_Esp, r2_Esp)


## ------------------------------------------------------------------------------------------- ## 


library(terra)


r1 <- rast("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_1985_int.tif")

r2 <- rast("E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_2019_int.tif")

esp <- vect("C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA/DATA_/VECTOR/Area_estudo/espinhaco_total_Simplf_WGS84.shp")


r1_Esp <- mask(r1, esp)

r2_Esp <- mask(r2, esp)



changeMatrixTerra <- crosstab(r1, r2)

changeMatrixTerraEspinhaco <- crosstab(c(r1_Esp, r2_Esp))

saveRDS(changeMatrixTerraEspinhaco,"C:/Users/JG/Desktop/changeMatrixTerraEspinhaco-v1.rds")
write.csv(changeMatrixTerraEspinhaco,"C:/Users/JG/Desktop/changeMatrixTerraEspinhaco-v1.csv")

