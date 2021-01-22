

if("jsonlite" %in% rownames(installed.packages()) == FALSE) {
  install.packages("jsonlite")
} 

library(jsonlite)

spData <- jsonlite::fromJSON("http://apiv3.iucnredlist.org/api/v3/species/chioglossa%20lusitanica?token=dd3def5af389de71d6b6543e2f67a9ab84b0f4b8c4dc019ef70adbda12452a39")

# Categoria IUCN
spData$result$category

# Tendência populacional
spData$result$population_trend
