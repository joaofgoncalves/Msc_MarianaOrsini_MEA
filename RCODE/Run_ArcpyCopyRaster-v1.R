
for(yr in 1986:2019){
  
  cat('arcpy.CopyRaster_management(
  in_raster=\"E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019/reclass_mapbiomas_bbox_',yr,'.tif\", 
  out_rasterdataset=\"E:/Projects/MSc_MarianaOrsini/DATA/RASTER/LULC_1985_2019_int/reclass_mapbiomas_bbox_',yr,'_int.tif", 
  config_keyword=\"\", 
  background_value=\"\", 
  nodata_value=\"\", 
  onebit_to_eightbit=\"NONE\", 
  colormap_to_RGB=\"NONE\", 
  pixel_type=\"8_BIT_UNSIGNED\", 
  scale_pixel_value=\"NONE\", 
  RGB_to_Colormap=\"NONE\", 
  format=\"TIFF\", 
  transform=\"NONE\")', '\n\n', sep="")
  
}
