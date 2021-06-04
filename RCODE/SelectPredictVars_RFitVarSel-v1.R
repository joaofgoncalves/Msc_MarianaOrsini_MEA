

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


## -------------------------------------------------------------------------------------------- ##


setwd("C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA")

source("./RCODE/_LibMainFunctions-v1.R")


## -------------------------------------------------------------------------------------------- ##


# Table with groups of variables: STATIC (like climate) and dynamic 
# (annual land use and satellite indices)
varGroups = read_csv("./DATA_/varNames-v1.csv")

# Base directory holding the train data/matrix
baseDir <- "./DATA_/TABLES/_TRAIN_DATASETS/BIRDS_TODO_AGAIN/"

# The actual column range of predictor variables to use that is in the TRAIN DATASETS files
varRange <- 8:52


## ------------------------------------------------------------------------------------------- ##


# List all train data and matrices files

flTrain <- list.files(baseDir, pattern="_TrainData_", full.names = TRUE)
flMat   <- list.files(baseDir, pattern="_TrainMatrix_", full.names = TRUE)

flTrain <- flTrain[grepl(".rds$", flTrain)]
flMat   <- flMat[grepl(".rds$", flMat)]


# Run the variable selection algorith by species
#
#

for(i in 1:length(flTrain)){
  
  trainData <- read_rds(flTrain[i])
  trainData <- trainData %>% select(-WETL)
  trainMatrix <- read_rds(flMat[i])
  predVars    <- colnames(trainData)[varRange]
  
  spName <- trainData[1,"SpeciesName"]
  
  cat(green("\n\n*****************************************************\n"))
  cat(green("*** Processing species name:", spName, "\n"))
  cat(green("*****************************************************\n\n"))
  
  selVars <- multiRoundIterVarSel(trainData, 
                                  trainMatrix, 
                                  predVars, 
                                  ntree        = 1000,
                                  nmax         = "10-rule", 
                                  tenRuleNmax  = 10,
                                  thresh       = 0.75, 
                                  aggFun       = median,
                                  method       = "spearman",
                                  varGroups, 
                                  doNZV        = TRUE, 
                                  doHighCor    = TRUE, 
                                  highCorVal   = 0.95,
                                  pbar         = TRUE, 
                                  verbose      = TRUE)
  
  print(selVars)
  outFn <- paste(baseDir, gsub("\\ +","_",spName),"_SelectedVars_v1.rds",sep="")
  write_rds(selVars,file = outFn)
}

