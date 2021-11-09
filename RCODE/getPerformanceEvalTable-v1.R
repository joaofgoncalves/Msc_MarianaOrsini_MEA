





library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)

#rm(list=ls())

spGroups <- list.dirs("./OUT/MODS/R02", recursive = FALSE)

selVarsDF <- readxl::read_excel("./DATA_/TABLES/LULC_Habitats_Species-simple-v2.xlsx")

## ---------------------------------------------------------------------------------- ##



k <- 0
for(spGroup in spGroups){
  
  spNames <- list.dirs(spGroup, recursive = FALSE) 
  
  spGroup_ <- basename(spGroup)
  
  
  for(spName in spNames){
    
    k <- k + 1
    spName_ <- basename(spName)
    
    spNameFull <- selVarsDF$Species[selVarsDF$SpeciesComp %in% spName_]
    
    
    spFileList <- list.files(spName, pattern="_EnsMod_evalDF_AllMetrics.csv$", full.names = TRUE)
    
    spNameFull <- selVarsDF$Species[selVarsDF$SpeciesComp %in% spName_]
    
    tmp <- read_csv(spFileList)[,1:5] %>% 
      `colnames<-`(c("Stat","StatValue","Threshold","Sensitivity","Specificity")) %>% 
    mutate(SpeciesName=spNameFull, SpeciesGroup=spGroup_) %>% 
      select(7,6,1:5)
    
    if(k==1){
      modEvalStats <- tmp
    }else{
      modEvalStats <- bind_rows(modEvalStats, tmp)
    }
    
  }
  
}

modEvalStats <- modEvalStats %>% filter(Stat != "KAPPA") %>% arrange(Stat,SpeciesGroup,SpeciesName)
modEvalStats[modEvalStats$Stat=="ROC","Stat"] <- "AUC"


write_csv(modEvalStats,"./OUT/MODS/R02/modPerformanceStats_AUC_TSS-R02-v1.csv")




