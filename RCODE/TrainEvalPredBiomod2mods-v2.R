

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
#outDir <- "E:/Projects/MSc_MarianaOrsini/OUT/MODS/R01"
outDir <- "./OUT/MODS/R02"


setwd(projectDir)

# Read ancillary functions
source("./RCODE/_LibMainFunctions-v1.R")

## -------------------------------------------------------------------------------------- ##

inputPath <- "/DATA_/TABLES/_TRAIN_DATASETS/v2/AMPHIBIANS/"


# Read 1-Km grid data in vector and raster
grid1k <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid.shp")
rstGrid1K <- raster("./DATA_/RASTER/Grid1K_IDs_WGS84UTM23S_v1.tif")

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
projNames <- paste("yr",yrs,sep="")

# List files from prediction for prediction datasets
flPredData <- list.files(paste(projectDir,"/DATA_/TABLES/_PRED_DATASETS/",sep=""), pattern = ".rds$", 
                         full.names = TRUE)
flPredData <- flPredData[grepl(paste(yrs, collapse="|"), flPredData)]


## -------------------------------------------------------------------------------------- ##

# Set the working directory for storing the models per species
setwd(outDir)

## -------------------------------------------------------------------------------------- ##


for(i in 1:length(flTrain)){
  
  # Load train data for each species from previously generated datasets
  trainData   <- read_rds(flTrain[i])
  trainData   <- trainData %>% select(-WETL)
  
  # Load the training matrix identifying presences and PA sets
  trainMatrix <- read_rds(flMat[i])
  
  # Load the predictive variables from RF iterative selection
  predVars    <- read_rds(flSelVars[i])
  
  
  
  # Create a 0/1 response variable
  respVar     <- unlist(trainData[,"pa"])
  
  # Create the species name
  spNameInit  <- trainData[1,"SpeciesName"]
  spName_     <- gsub("\\ +","_",spNameInit)
  spName      <- gsub("\\ +","",spNameInit)
  
  
  cat(green("\n-> Processing species name:", spName, "\n\n"))
  
  # Format data from trainData according to biomod2 standards
  myBiomodData <- BIOMOD_FormatingData(resp.var       = respVar,
                                       expl.var       = trainData[,predVars],
                                       resp.name      = spName,
                                       PA.strategy    = 'user.defined',
                                       PA.table       = trainMatrix) # PA generation method
  
  # Define some modelling options for available algorithms
  myBiomodOptions <- BIOMOD_ModelingOptions(GAM = list(k = 4),
                                            MAXENT.Phillips = list(threshold=FALSE,
                                                                   hinge=FALSE,
                                                                   path_to_maxent.jar="C:/MyFiles/temp"), ## Change this too!!!!!!!
                                            GBM = list(n.trees = 2000))
  
  
  ## -------------------------------------------------------------------------------------- ##
  ## Calibrate models ----
  ## -------------------------------------------------------------------------------------- ##
  
  # Train models across PA sets and evaluation rounds
  #
  #
  myBiomodModelOut <- suppressMessages(BIOMOD_Modeling(
    data           = myBiomodData, # Input data
    models = c('GLM','GBM','GAM',
               'CTA','ANN','FDA',
               'MARS','RF','MAXENT.Phillips.2'), # Models to run
    models.options = myBiomodOptions,
    NbRunEval      = 10, # Number of Evaluation runs
    DataSplit      = 75,
    Prevalence     = 0.5, # Prevalence between 0 and 1
    VarImport      = 5, # Nr of rounds to evaluate variables
    models.eval.meth  = c('TSS','ROC','KAPPA'), # Evaluation metrics
    SaveObj           = TRUE,
    rescal.all.models = TRUE,
    do.full.models    = TRUE)) # Model with all data?
  
  # Get model evaluation values
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # Print ROC scores
  print(myBiomodModelEval["ROC","Testing.data",,,])
  print(myBiomodModelEval["TSS","Testing.data",,,])
  
  # Get boxplot stats
  print(fivenum(as.numeric(myBiomodModelEval["ROC","Testing.data",,,])))
  print(fivenum(as.numeric(myBiomodModelEval["TSS","Testing.data",,,])))
  
  # Save evaluation metrics from the arrays
  evalDF.ROC <- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
  evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
  evalDF.KAPPA <- as.data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
  
  write.csv(evalDF.ROC, file = paste(getwd(),"/",spName,"/",spName,"_evalDF_ROC.csv",sep=""))
  write.csv(evalDF.TSS, file = paste(getwd(),"/",spName,"/",spName,"_evalDF_TSS.csv",sep=""))
  write.csv(evalDF.KAPPA, file = paste(getwd(),"/",spName,"/",spName,"_evalDF_KAPPA.csv",sep=""))
  
  
  ## -------------------------------------------------------------------------------------- ##
  ## Perform ensemble modelling ----
  ## -------------------------------------------------------------------------------------- ##
  
  
  # Select best models in a two-step fashion
  # First select the best six algorithms (out of nine) and then choose the top 10% models 
  # for each technique
  selMods <- twoStepBestModelSelection(myBiomodModelOut, 
                                       evalMetric = "TSS", 
                                       nrBestAlgos = 6, 
                                       bestAlgoFun = stats::median, 
                                       topFraction = 0.1)

    
  # Perform the ensemble of best models previously selected
  myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                        chosen.models = selMods,
                                        em.by = 'all',
                                        prob.mean = TRUE,
                                        prob.mean.weight.decay = 'proportional')
  
  # Get evaluation scores for the Ensemble Modeling stage
  emEvalDF <- as.data.frame(get_evaluations(myBiomodEM))
  write.csv(emEvalDF, file = paste(getwd(),"/",spName,"/",spName,
                                   "_EnsMod_evalDF_AllMetrics.csv",sep=""))
  
  
  ## -------------------------------------------------------------------------------------- ##
  ## Obtain spatiotemporal projections for each year ----
  ## -------------------------------------------------------------------------------------- ##
  
  # Models to consider in the ensemble and projection
  modelsToUse <- get_kept_models(myBiomodEM, 1)
  
  # k=1
  # projName <- projNames[k]
  
  k <- 0
  for(projName in projNames){
    
    k <- k + 1
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
  
  T# Save the accumulated projections data frame
  write_rds(allProjs,paste("./",spName,"/",spName,"_allProjections-v1.rds",sep=""))
}




