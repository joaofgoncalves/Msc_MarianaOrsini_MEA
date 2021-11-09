

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(exactextractr)
library(stringi)
library(stringr)
library(rlang)
library(crayon)
library(randomForest)
library(biomod2)
library(biomod2plus)
library(readxl)

projectDir <- "C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA"

inputPaths <- c("./DATA_/TABLES/_TRAIN_DATASETS/v2/AMPHIBIANS",
                "./DATA_/TABLES/_TRAIN_DATASETS/v2/BIRDS",
                "./DATA_/TABLES/_TRAIN_DATASETS/v2/REPTILES")



## -------------------------------------------------------------------------------------- ##


setwd(projectDir)

# Read ancillary functions
source("./RCODE/_LibMainFunctions-v1.R")


## -------------------------------------------------------------------------------------- ##

k<-0
for(inputPath in inputPaths){
  
  spGroup <- basename(inputPath)
  
  # List files containing selected variables from previous steps
  flSelVars <- list.files(paste(projectDir,inputPath,sep=""), pattern="_SelectedVars_", 
                          full.names = TRUE)
  flSelVars <- flSelVars[grepl(".rds$", flSelVars)]
  
  
  selVarsDF <- read_excel("./DATA_/TABLES/LULC_Habitats_Species-simple-v2.xlsx")
  
  
  ## -------------------------------------------------------------------------------------- ##
  
  for(i in 1:length(flSelVars)){
    k<-k+1
    # Create the species name
    spNameInit  <-     paste(unlist(strsplit(basename(flSelVars[i]),"_"))[1:2], collapse=" ")
    spName_     <- gsub("\\ +","_",spNameInit)
    spName      <- gsub("\\ +","",spNameInit)
    
    
    
    # Load the predictive variables from RF iterative selection
    predVars    <- read_rds(flSelVars[i])
    
    spIdx <- trimws(selVarsDF$Species, which = "both") %in% spNameInit
    
    selVarsExp <- selVarsDF[spIdx,c("H1","H2","H2.1","H3")] %>% as.character() %>% na.omit()
    
    predVars <- unique(c(predVars, selVarsExp))
    
    tmp <- data.frame(SpeciesGroup = spGroup, 
                      SpeciesName  = spNameInit,
                      nVars        = length(predVars),
                      PredVars     = paste(predVars,collapse=", "))
    
    if(k==1){
      preVarsUsed <- tmp
    }else{
      preVarsUsed <- bind_rows(preVarsUsed,tmp)
    }
    
  }
}


write_csv(preVarsUsed,"./OUT/MODS/R02/preVarsUsed-R02-v1.csv")




