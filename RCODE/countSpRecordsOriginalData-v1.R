

library(readxl)
library(dplyr)


spData_anf <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Amphibians/endemic_amphibia_Espinhaco_v2-1.xlsx")

spData_bir <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Birds/DH_Highland_Birds-v3.xlsx")

spData_rep <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Reptiles/Repteis-Espinhaco_Database_2021-vii-15-Henrique_MO_JG_v3.xlsx")


spCounts <- rbind(
  
  spData_anf %>% 
  mutate(spGroup = "Amphibians") %>% 
  group_by(spGroup, Species) %>% 
  summarise(spCount = n()) %>% 
  as.data.frame() %>% 
  arrange(Species),
  
  spData_bir %>% 
    rename(Species = SpeciesName) %>% 
    mutate(spGroup = "Birds") %>% 
    group_by(spGroup, Species) %>% 
    summarise(spCount = n()) %>% 
    as.data.frame() %>% 
    arrange(Species),
  
  spData_rep %>% 
    mutate(spGroup = "Reptiles") %>% 
    group_by(spGroup, Species) %>% 
    summarise(spCount = n()) %>% 
    as.data.frame() %>% 
    arrange(Species)
  
)
  
write.csv(spCounts,"./OUT/SpeciesCounts-OriginalRecords-v1.csv", row.names = FALSE)
