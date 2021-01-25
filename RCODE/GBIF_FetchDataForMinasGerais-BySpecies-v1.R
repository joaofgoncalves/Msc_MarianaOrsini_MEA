


st_envelope <- function(x) st_as_sfc(st_bbox(x))

shpMinasGerais <- read_sf("./DATA_/VECTOR/gadm36_BRA_shp/gadm36_BRA_1.shp") %>% 
  filter(NAME_1 == "Minas Gerais") %>% 
  st_envelope()

wktMG <- st_geometry(shpMinasGerais[1,]) %>% st_as_text()

library(rgbif)
library(sf)
library(dplyr)

# Area de estudo de Minas Gerais
wktMG <- "POLYGON ((
-51.04588 -22.92275, 
-39.85676 -22.92275, 
-39.85676 -14.23343, 
-51.04588 -14.23343, 
-51.04588 -22.92275))"

dtSpecies <- occ_data(scientificName = "Platalea ajaja",
                      limit = 100000, geometry = wktMG)
