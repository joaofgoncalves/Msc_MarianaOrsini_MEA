

library(readxl)
library(readr)
library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(rgbif)
library(lubridate)
library(ggmap)
library(OpenStreetMap)
library(ggplot2)


#spList <- read_excel("DATA_/SP_DATA/Potenciais_Sp_Alvo_rev.xlsx")
spList <- read_excel("DATA_/SP_DATA/Potenciais_Sp_Alvo_rev.xlsx", sheet = 1)

bbox <- read_sf("./DATA_/TEMP/SerraEspinhaco_bbox_50Km.shp")

buff0 <- read_sf("./DATA_/TEMP/SerraEspinhaco_buff50k_bbox.shp") %>% 
  st_transform(crs="EPSG:4326")
#plot(buff)


buff <- read_sf("./DATA_/VECTOR/SerraEspinhaco/SerraEspinhaco_Dissolve_v3.shp") %>% 
  st_transform(crs="EPSG:4326")  
  #st_buffer(1E-3)
#plot(buff)

coordsDF <- st_coordinates(buff) %>% 
  as.data.frame
coordsBuff0DF <- st_coordinates(buff0) %>% 
  as.data.frame


wktSE <- st_geometry(bbox) %>% st_as_text()
wktCentroid <- st_geometry(bbox) %>% st_centroid() %>% st_as_text()


mapOSM <- openmap(upperLeft  = c(-9.307559, -45.305857),
                  lowerRight = c(-21.310615, -39.094775), type="stamen-terrain", zoom=7)



map.latlon <- openproj(mapOSM, projection = "EPSG:4326")




i <- 0
pb <- txtProgressBar(min=1, max=length(spList$Species_name), style=3)

dtSpecies<-data.frame()

for(spName in spList$Species_name){
  
  i <- i + 1
  
  dtSpeciesTmp <- occ_data(scientificName = spName,
                        limit = 1000000, geometry = wktSE)$data 

  if(!is.null(dtSpeciesTmp)){

    
    ptsDF <- dtSpeciesTmp[,c("decimalLongitude", "decimalLatitude")]
    pts <- st_multipoint(as.matrix(ptsDF))


    subTt <- paste(dtSpeciesTmp[1,c("order","family","genus")],collapse = ", ")

    mapIt <- OpenStreetMap::autoplot.OpenStreetMap(map.latlon, plot=FALSE) +
      #geom_sf(data=buff, fill=NA) +
      geom_polygon(data = coordsDF, mapping = aes(x=X, y=Y), fill=NA, color="black") +
      geom_polygon(data = coordsBuff0DF, mapping = aes(x=X, y=Y), fill=NA, color="black") +
      geom_point(data=ptsDF, mapping=aes(x=decimalLongitude, y=decimalLatitude), fill=NA,
                 color = "red", size=3) +
      coord_sf(crs="EPSG:4326") +
      labs(title = paste("[",i,"] ",spName," (n=",nrow(dtSpeciesTmp),")",sep=""), subtitle = subTt) +
      xlab("Longitude") +
      ylab("Latitude")

    ggsave(filename=paste("./OUT/spRecordsGBIF_maps/",i,"_",spName,".png",sep=""),
           mapIt,width = 7, height = 12)
    
    
    if(i==1){
      dtSpecies <- dtSpeciesTmp
    }else{
      dtSpecies <- bind_rows(dtSpecies, dtSpeciesTmp)
    }
  }
  setTxtProgressBar(pb, i)
}


write.csv(dtSpecies,"./DATA_/TABLES/dtSpecies-v1.csv", row.names = FALSE)


dtSpeciesSel$species %>% unique

dtSpeciesSel <- 
dtSpecies %>%
  dplyr::select(key,
                scientificName,
                phylum,
                order,
                family,
                genus,
                species,
                eventDate,
                dateIdentified,
                decimalLatitude,
                decimalLongitude,
                coordinateUncertaintyInMeters)


write.csv(dtSpeciesSel,"./DATA_/TABLES/dtSpeciesSel-v1.csv", row.names = FALSE)



dtSpeciesSel %>% 
  mutate(eventDate1 = as.Date(eventDate)) %>% 
  filter(year(eventDate1) >= 1980) %>% 
  filter(coordinateUncertaintyInMeters < 2500 | is.na(coordinateUncertaintyInMeters)) %>% 
  group_by(species) %>% 
  summarize(nOcc=n()) %>% 
  filter(nOcc > 20) %>% 
  left_join(spList, by=c("species"="Species_name")) %>%
  arrange(Group) %>% 
  View

