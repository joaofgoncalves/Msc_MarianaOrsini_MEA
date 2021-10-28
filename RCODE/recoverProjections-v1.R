
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

outDir <- "D:/Projects/MSc_MarianaOrsini/OUT/MODS/R02/BIRDS"

inputPath <- "./DATA_/TABLES/_TRAIN_DATASETS/v2/BIRDS/"



## -------------------------------------------------------------------------------------- ##


setwd(projectDir)

# Read ancillary functions
source("./RCODE/_LibMainFunctions-v1.R")


## -------------------------------------------------------------------------------------- ##

clNamesDF <- read_excel("./DATA_/TABLES/MapBiomasNewClasses-v1.xlsx")
clNames <- clNamesDF %>% pull(ClassCode)



## -------------------------------------------------------------------------------------- ##

# Read 1-Km grid data in vector and raster
grid1k <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid.shp")
rstGrid1K <- raster("./DATA_/RASTER/Grid1K_IDs_WGS84UTM23S_v1.tif")


## -------------------------------------------------------------------------------------- ##



# List datasets used for training
flTrain <- list.files(paste(projectDir,inputPath,sep=""), pattern="_TrainData_", 
                      full.names = TRUE)
flTrain <- flTrain[grepl(".rds$", flTrain)]

# List matrices used for train with info on Pseuso-absence sets
flMat <- list.files(paste(projectDir,inputPath,sep=""), pattern="_TrainMatrix_", 
                    full.names = TRUE)
flMat <- flMat[grepl(".rds$", flMat)]

# List files containing selected variables from previous steps
flSelVars <- list.files(paste(projectDir,inputPath,sep=""), pattern="_SelectedVars_", 
                        full.names = TRUE)
flSelVars <- flSelVars[grepl(".rds$", flSelVars)]


selVarsDF <- read_excel("./DATA_/TABLES/LULC_Habitats_Species-simple-v2.xlsx")

## -------------------------------------------------------------------------------------- ##



# Generate year sequence
yrs <- seq(1985,2019, by=2)

# Projectio names
projNamesInit <- paste("yr",yrs,sep="")

# List files from prediction for prediction datasets
flPredData <- list.files(paste(projectDir,"/DATA_/TABLES/_PRED_DATASETS/",sep=""), pattern = ".rds$", 
                         full.names = TRUE)
flPredData <- flPredData[grepl(paste(yrs, collapse="|"), flPredData)]


spName <- "Cinclodesespinhacensis"
spNameSep <- "Cinclodes_espinhacensis"
spNameInit <- "Cinclodes espinhacensis"

load("./OUT/MODS/R02/BIRDS/Cinclodesespinhacensis/Cinclodesespinhacensis.1635301033.models.out")
load("./OUT/MODS/R02/BIRDS/Cinclodesespinhacensis/Cinclodesespinhacensis.1635301033ensemble.models.out")

myBiomodModelOut <- get(ls()[grepl(".models.out",ls())][1])
myBiomodEM <- get(ls()[grepl(".models.out",ls())][2])


projNames <- projNamesInit[9:18]


## -------------------------------------------------------------------------------------- ##

# Set the working directory for storing the models per species
setwd(outDir)

## -------------------------------------------------------------------------------------- ##

idx <- 1:length(flSelVars)
idx <- idx[grepl(spNameSep,flSelVars)]

predVars    <- read_rds(flSelVars[idx])
spIdx <- trimws(selVarsDF$Species, which = "both") %in% spNameInit

if(!any(spIdx)){
  cat(red("\n\n-------------------------------------------------------\n"))
  cat(red("Did not found species name:\n", spNameInit,"\n\n"))
}

selVarsExp <- selVarsDF[spIdx,c("H1","H2","H2.1","H3")] %>% as.character() %>% na.omit()

predVars <- unique(c(predVars, selVarsExp))



## -------------------------------------------------------------------------------------- ##
## Obtain spatiotemporal projections for each year ----
## -------------------------------------------------------------------------------------- ##

# Models to consider in the ensemble and projection
modelsToUse <- get_kept_models(myBiomodEM, 1)

# k=1
# projName <- projNames[k]

k <- 0
for(projName in projNamesInit){
  
  k <- k + 1
  
  if(!(projName %in% projNames)){
    cat("Skipping projection name:", projName,"\n\n")
    next
  }
  
  predData <- read_rds(flPredData[k])
  
  # Obtain spName spatiotemporal projections for all the best selected models
  myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                    new.env         = predData[,predVars],
                                    proj.name       = projName,
                                    selected.models = modelsToUse,
                                    filtered.meth   = NULL,
                                    binary.meth     = NULL,
                                    compress        = 'gzip',
                                    clamping.mask   = TRUE,
                                    output.format   = '.RData',
                                    do.stack        = TRUE,
                                    keep.in.memory  = FALSE,
                                    on_0_1000       = TRUE)
  
  # Perform the ensembling of projections
  myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                           binary.meth       = c('TSS','ROC','KAPPA'),
                                           EM.output         = myBiomodEM,
                                           output.format     = '.RData')
  
  # Rename the columns of the projections datasets
  cnames1 <- paste(projName, c('TSSbin','ROCbin','KAPPAbin'),sep="_")
  ensProjMatrix <- getEnsembleBinProjs(spName, projName)
  colnames(ensProjMatrix) <- cnames1
  
  cnames2 <- paste(projName,"HS_ens",sep="_")
  ensProjHsMatrix <- getEnsembleHsProjs(spName, projName)[,1,drop=FALSE]
  colnames(ensProjHsMatrix) <- cnames2
  
  # Create a table with projections for the entire reference grid       
  ensProjMatrix <- bind_cols(predData   %>% dplyr::select(ID),
                             ensProjHsMatrix %>% dplyr::select(1),
                             ensProjMatrix) %>% suppressMessages()
  
  # Assign spatial to the projections grid
  tmpGrid <- grid1k %>% left_join(ensProjMatrix, by="ID")
  
  # Convert the binary and habitat suitability projections data from vector to raster
  tmpRst_Bin <- fasterize::fasterize(tmpGrid, raster = rstGrid1K, field = cnames1[1])
  tmpRst_Hsb <- fasterize::fasterize(tmpGrid, raster = rstGrid1K, field = cnames2[1])
  
  # Write raster data for the binary and HS projections by year
  fnBin <- paste("./",spName,"/",spName,"_TSSbin_",projName,".tif",sep="")
  fnHsb <- paste("./",spName,"/",spName,"_hsb_",projName,".tif",sep="")
  writeRaster(tmpRst_Bin, filename = fnBin, datatype="INT1U")
  writeRaster(tmpRst_Hsb, filename = fnHsb, datatype="INT2U")
  
  # Accumulate all annual binary and HS projections for the reference grid in one large data frame
  if(k==1){
    allProjs <- ensProjMatrix
  }else{
    allProjs <- allProjs %>% left_join(ensProjMatrix, by="ID")
  }
}

# Save the accumulated projections data frame
write_rds(allProjs,paste("./",spName,"/",spName,"_allProjections-v1.rds",sep=""))






