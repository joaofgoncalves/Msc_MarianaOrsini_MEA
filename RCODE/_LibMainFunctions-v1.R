
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



## Get functions for returning basic inputs -----------------------------------

#' Get biomod2 binary projections
#' 
#' Reads the RData file containing ensemble projections with a file
#' suffix similar to _ensemble_TSSbin.RData
#' 
#' @param spName The species name
#' @param projName Projection name
#' 
#' @return A data.frame with binary ensemble projections. 
#' Each column corresponds to a different combination of ensembling 
#' method / performance metric
#' 

getEnsembleBinProjs <- function(spName, projName){
  
  outFn <- paste("./",spName,"/proj_",projName,"/","proj_",projName,"_",
                 spName,"_ensemble_TSSbin.RData",sep="")
  load(outFn)
  return(as.data.frame(get(paste("proj",projName,spName,"ensemble_TSSbin",sep="_"))))
  
}


#' Get biomod2 habitat suitability (continuous) projections
#' 
#' Reads the RData file containing ensemble projections with a file
#' suffix equal to _ensemble.RData
#' 
#' @param spName The species name
#' @param projName Projection name
#' 
#' @return A data.frame with habitat seuitability ensemble projections. 
#' Each column corresponds to a different combination of ensembling 
#' method / performance metric
#' 

getEnsembleHsProjs <- function(spName, projName){
  
  outFn <- paste("./",spName,"/proj_",projName,"/","proj_",projName,"_",
                 spName,"_ensemble.RData",sep="")
  load(outFn)
  #return(as.data.frame(get(paste("proj",projName,spName,"ensemble",sep="_"))))
  return(as.data.frame(ef.out))
}

#' Get unique years in a dataset
#' 
#' An ancillary function to retrieve a vector with all different/unique 
#' years in a species dataset
#'  
#' @param x A data.frame with a column named "Year" prepared by 
#' \code{prepSpData} function
#' 
#' @return A sorted vector with years
#' 

getYears <- function(x){
  return(x %>% pull(Year) %>% unique %>% sort)
} 


#' Get unique species names in a dataset
#' 
#' An ancillary function to retrieve a vector with all different/unique 
#' species names in an occurrence/presence-only dataset
#'  
#' @param x A data.frame with a column named "SpeciesName" prepared by 
#' \code{prepSpData} function
#' 
#' @return A sorted vector with species names
#' 

getSpeciesNames <- function(x){
  return(x %>% pull(SpeciesName) %>% unique %>% sort)
}

#' Get unique ID's in a dataset
#' 
#' An ancillary function to retrieve a vector with all different/unique 
#' IDs (integers) in an occurrence/presence-only dataset
#'  
#' @param x A data.frame with a column named "ID" prepared by 
#' \code{prepSpData} function
#' 
#' @return A sorted vector with IDs
#' 

getIDs <- function(x){
  return(x %>% pull(ID))
} 

#' Get unique ID's per year in a dataset
#' 
#' An ancillary function to retrieve a vector with all different/unique 
#' IDs (integers) in an occurrence/presence-only dataset
#'  
#' @param x A data.frame with a column named "ID" prepared by 
#' \code{prepSpData} function
#' 
#' @param yr Target year (integer)
#' 
#' @return A sorted vector with IDs
#' 

getYearIDs <- function(x, yr){
  return(x %>% filter(Year == yr) %>% pull(ID))
}


#' Read data files by year
#' 
#' From a list of datasets with the year in filename it picks a 
#' specific dataset for a target year and reads it
#'  
#' @param fileList A vector with file paths to annual datasets
#' 
#' @param year Target year (integer)
#' 
#' @return A data.frame/tibble object
#' 


## Functions for reading, processing and preparing data for model development -----


readDatasetByYear <- function(fileList, year){
  fpath <- fileList[grepl(year,fileList)]
  return(read_rds(fpath))
}

#' Merge data frames in a list
#' 
#' A simple function
#'  
#' @param x A list with one data.frame per slot
#' 
#' @return A data frame merging all partial data
#' 

mergeDFsInList <- function(x){
  for(i in 1:length(x)){
    if(i==1){
      outDF <- x[[i]]
    }else{
      outDF <- bind_rows(outDF, x[[i]])
    }
  }
  return(outDF)
}


#' Prepare species data of occurrences
#' 
#' Applies a filter a renaming to species occurrence data 
#' preparing it for use in models
#'  
#' @param spData A data.frame with species data
#' @param spNameCol Column name with species names
#' @param lonCol Column name with longitude data
#' @param latCol Column name with latitude data
#' @param yearCol Column name with year data
#' @param yrStart Start year of the target period
#' @param yrEnd End year of the target period
#' @param filterByMinObs Filter species by a minimum number of 
#' observations? (default: TRUE)
#' @param nmin Minimum number of observations
#' @param asSF Return object as a simple features data object 
#' (sf package, fefault: TRUE)
#' 
#' @return A data.frame/tibble object with prepared species data
#' 

prepSpData <- function(spData, 
                       spNameCol      = "Species", 
                       lonCol         = "Lon", 
                       latCol         = "Lat", 
                       yearCol        = "Year",
                       yrStart        = 1985, 
                       yrEnd          = 2019, 
                       filterByMinObs = TRUE,
                       nmin           = 30,
                       asSF           = TRUE){
  
  outDF <- spData %>% 
    select(!!spNameCol, !!lonCol, !!latCol, !!yearCol) %>% 
    filter(.[[yearCol]] >= yrStart & .[[yearCol]] <= yrEnd) %>% 
    arrange(across(.cols = all_of(c(spNameCol, yearCol)))) %>% 
    `colnames<-`(c("SpeciesName","Lon","Lat","Year"))
  
  if(filterByMinObs){
    outDF <- outDF %>% 
    group_by(SpeciesName) %>% 
      filter(n() >= nmin) %>% 
      ungroup() %>% 
      arrange(SpeciesName, Year)
    
  }
  
  if(asSF){
    
    spPtDF <- SpatialPointsDataFrame(outDF[,c("Lon","Lat")], 
                                     data=outDF, 
                                     proj4string = CRS("EPSG:4326"))
    
    return(st_as_sf(spPtDF))
    
  }else{
    return(outDF)
  }
}

#' Get grid identifiers intersecting species occurrence points
#' 
#' Given a species occurrence dataset the function applies a step to 
#' retrieve the spatial intersection between species point records and a 
#' reference grid in vector format. It then performs duplicate removal and 
#' filters species whose minimum number of records is below a user-defined 
#' threshold.
#' 
#' @param sfSpeciesDF A sf object generated by \code{prepSpData}
#' @param sfGrid A sf object with the reference grid
#' @param removeDups Remove duplicates? (default: TRUE)
#' @param filterByMinObs Filter by a minimum number of observations? (default: TRUE)
#' @param nmin Mimimum number of gridded observations (not raw points) threshold
#' @param getCounts Calculate species counts (i.e. gridded observations)
#' 
#' @return If getCounts=FALSE returns a sf object with grid data including the ID 
#' and the geometry. If TRUE then returns a list with two elements one for the species 
#' data ("spData") and the other one with a species count table ("spCounts")
#' 
#' @note A  gridded observation is considered duplicated if it has the same year.
#' Multiple observations for the same quadrat/grid element are allowed if they 
#' are recorded for different years within the target period.
#' 

prepSpDataWithGridID <- function(sfSpeciesDF, 
                       sfGrid, 
                       removeDups     = TRUE,
                       filterByMinObs = TRUE,
                       nmin           = 30,
                       getCounts      = FALSE){
  
  outDF <- st_intersection(sfSpeciesDF, sfGrid) %>% 
    suppressMessages() %>% 
    suppressWarnings() %>% 
    na.omit() %>% 
    arrange(SpeciesName, ID, Year)
  
  if(removeDups){
    
    # Remove records by repeated 1km grid element
    # for the same year
    # Different years for same species / ID are kept!!!
    
    outDF <- outDF %>% 
      mutate(cd = paste(gsub("\\ +","_",SpeciesName),ID,Year,sep="_")) %>% 
      filter(!duplicated(cd))
      #distinct(SpeciesName, ID, Year, 
      #                         .keep_all = TRUE) # DOES NOT WORK
  }
  
  if(filterByMinObs){
    # Filter by small datasets
    outDF <- outDF %>% 
      group_by(SpeciesName) %>% 
      filter(n() >= nmin) %>% 
      ungroup() %>% 
      arrange(SpeciesName, Year)
  }

  if(getCounts){
    aggDF <- outDF %>%
      #st_drop_geometry() %>%
      group_by(SpeciesName) %>%
      summarize(nGridCount = n()) %>%
      arrange(desc(nGridCount))
    
    return(list(spData   = outDF,
                spCounts = aggDF))
    
  }else{
    return(outDF)
  }
}

#' Get environmental data by year
#' 
#' An ancillary function that retrieves the environmental data for each species 
#' records on an annual basis. This means that data is fetched by each year block to 
#' match observations temporally with 'dynamic' environmental data. 
#' The data is retrieved by joining species records and environmental data by 
#' the identifier code of each quadrat/grid element. The function works for multiple 
#' species data.
#' 
#' @param dfPathList A list of file paths containing the full annual 
#' datasets used for prediction. Each data frame must include all predictive variables 
#' mapped for every grid element
#' 
#' @param spDataGrid A sp data object generated by \code{prepSpDataWithGridID}
#' 
#' @param progressBar Print a progress bar? (default: TRUE)
#' 
#' @return a data frame with env data matching gridded locations of species records

getEnvDataByYear <- function(dfPathList,
                             spDataGrid,
                             progressBar = TRUE){

  if("spData" %in% names(spDataGrid)){
    spDataGrid <- spDataGrid[["spData"]]
  }
  
  if(inherits(spDataGrid,"sf")){
    spDataGrid <- spDataGrid %>% st_drop_geometry()
  }

  spDataGrid <- data.frame(pa = 1, as.data.frame(spDataGrid))
  
  yrs <- getYears(spDataGrid)
  
  if(progressBar){
    pb <- txtProgressBar(1,length(yrs),style = 3)
    message(green("\nReading presence data by year:\n"))
  }
  
  for(i in 1:length(yrs)){
    
    spDataGridTmp <- spDataGrid %>% filter(Year == yrs[i])
    envDataYr <- readDatasetByYear(dfPathList, yrs[i])
    
    tmpDF <- spDataGridTmp %>% 
      left_join(envDataYr, by = "ID")
    
    if(i==1){
      annualDF <- tmpDF
    }else{
      annualDF <- bind_rows(annualDF, tmpDF)
    }
    
    if(progressBar) setTxtProgressBar(pb,i)
  }

  return(annualDF %>% arrange(SpeciesName, Year, ID))
}

#' Generate random pseudo-absences
#' 
#' A function used to generate random pseudo-absences in an annual basis, 
#' i.e., matching the year of presence records and at the same quantity. 
#' 
#' @param spName The species name
#' @param dfPathList A vector with file paths for the full prediction 
#' datasets with environmental data 
#' @param envDataOcc Data generated by the \code{getEnvDataByYear} function
#' @param nPAsets Number of pseudo-absences sets (default: 10)
#' @param nPAperSet Number of pseudo-absences per set (default: "equal"). 
#' If set to "equal", the number of PA's is set to the number of presences 
#' in each PA set
#' @param progressBar Print a progress bar? (default: TRUE)
#' 
#' @return A list object with each set of pseudo absences per slot
#' 

doRandomPA <- function(spName,
                        dfPathList,
                        envDataOcc,
                        nPAsets    = 10,
                        nPAperSet  = "equal",
                        progressBar = TRUE){
  
  envDataOccTargetSp <- 
    envDataOcc %>% filter(SpeciesName == spName)
  
  yrs    <- getYears(envDataOccTargetSp)
  ids    <- getIDs(envDataOccTargetSp)
  allIDs <- getIDs(envDataOcc)

  if(nPAperSet == "equal"){
    nPAperSet <- nrow(envDataOccTargetSp)
  }
  
  nPAperYear <- round(nPAperSet / length(yrs))
  
  PseudoAbsData <- list()
  
  if(progressBar){
    pb <- txtProgressBar(1,length(yrs)*nPAsets,style = 3)
    #message(green("\nGenerating pseudo-absences data by year:\n"))
  }
  
  
  z <- 0
  k <- 0
  for(yr in yrs){
    
    envDataYr <- readDatasetByYear(dfPathList, yr) %>% 
      filter(!(ID %in% allIDs)) %>% 
      na.omit() %>% 
      as_tibble()
    
    k <- k + 1 # Year counter
    
    for(i in 1:nPAsets){
      z <- z + 1
      tmpPAdf <- data.frame(pa = 0, 
                            sample_n(envDataYr, nPAperYear))
      
      if(k==1){
        PseudoAbsData[[i]] <- tmpPAdf
      }else{
        PseudoAbsData[[i]] <- bind_rows(PseudoAbsData[[i]], tmpPAdf)
      }
      
      if(progressBar) setTxtProgressBar(pb, z)
      
    }
    
  }

  attr(PseudoAbsData,"SpeciesName") <- spName
  attr(PseudoAbsData,"nPresences")  <- nrow(envDataOccTargetSp)
  attr(PseudoAbsData,"nPAsets")     <- nPAsets
  attr(PseudoAbsData,"nPAperSet")   <- nPAperSet
  attr(PseudoAbsData,"nPAperYear")  <- nPAperYear
  attr(PseudoAbsData,"yrs")         <- yrs
  attr(PseudoAbsData,"removedIds")  <- allIDs
  
  return(PseudoAbsData)
}

#' Create the train matrix
#' 
#' An ancillary function used to create a boolean train matrix suited for use 
#' in biomod2 package. This matrix has as many columns as PA sets. Value TRUE 
#' indicates the data to use for each set (both presences - at the top rows - 
#' and followed by PA's)
#' 
#'  @param PseudoAbsData An object created by the \code{doRandomPA} function
#'  
#'  @return A boolean train matrix suited for biomod2 

doTrainMatrix <- function(PseudoAbsData){
  
  nPresences <- attr(PseudoAbsData,"nPresences")
  nPAsets    <- attr(PseudoAbsData,"nPAsets") 
  nPAperYear <- attr(PseudoAbsData,"nPAperYear") 
  yrs        <- attr(PseudoAbsData,"yrs")
  
  ntotalPAs <- nPAperYear*length(yrs)*nPAsets
  ntoralRows <- ntotalPAs+nPresences
 
  trainMatrix <- matrix(FALSE, nrow = ntoralRows, ncol = nPAsets)
  trainMatrix[1:nPresences,] <- TRUE
  
  rowSt <- (nPresences + 1)
  
  for(i in 1:ncol(trainMatrix)){
    rowEnd <- (rowSt + (nPAperYear*length(yrs))) - 1
    trainMatrix[rowSt:rowEnd, i] <- TRUE
    rowSt <- rowEnd + 1
  }
  
  attr(trainMatrix,"nPresences") <- nPresences
  attr(trainMatrix,"nPAsets")    <- nPAsets
  attr(trainMatrix,"nPAperYear") <- nPAperYear
  attr(trainMatrix,"yrs")        <- yrs
  attr(trainMatrix,"ntotalPAs")  <- ntotalPAs
  attr(trainMatrix,"ntoralRows") <- ntoralRows
  
  #View(trainMatrix)
  return(trainMatrix)
  
}


#' Create train data from an initial dataset by species
#' 
#' This function is a wrapper that is used to:
#' 1. Select target columns from the initial data: latitude, longitude and year of the record;
#' 2. Rename columns
#' 3. Extract reference grid ID's (column name should be ID)
#' 4. Remove duplicates (if these are for the same grid ID and year, if records are for the 
#' same grid ID but different year these will be kept)
#' 5. Remove species if these have a small number of records
#' 6. Read environmental data for presence records using the year using as a basis 
#' and a set of annualized dataframes combining static and dynamic variables
#' 7. Generate pseudo-absence sets through random selection
#' 
#' 
#' @rdname prepSpData
#' @rdname getEnvDataByYear
#' @rdname getEnvDataByYear
#' @rdname doRandomPA
#' 
#' @param outDir An output directory to place output train files 
#' for each species
#' 
#' @return A data frame with output data. Also a set of files for each species with the train 
#' data and the train matrix to be used for modelling using biomod2's package 
#' 
#' @note A column named "pa" is added to the data showing presences (1's) and 
#' pseudo-absences (0's)
#' 

createTrainData <- function(spData, 
                            sfGrid,
                            dfPathList,
                            spName         = NULL,
                            spNameCol      = "Species", 
                            lonCol         = "Lon", 
                            latCol         = "Lat", 
                            yearCol        = "Year",
                            yrStart        = 1985, 
                            yrEnd          = 2019, 
                            filterByMinObs = TRUE,
                            nmin           = 30,
                            removeDups     = TRUE,
                            nPAsets        = 10,
                            nPAperSet      = "equal",
                            progressBar    = TRUE,
                            outDir         = getwd()
                            ){
  
  
  message(blue("\n\n***** PREPARING AND FILTERING DATA *****\n"))
  
  sfSpeciesDF <- prepSpData(
    spData         = spData,
    spNameCol      = spNameCol,
    lonCol         = lonCol,
    latCol         = latCol,
    yearCol        = yearCol,
    yrStart        = yrStart,
    yrEnd          = yrEnd,
    filterByMinObs = filterByMinObs,
    nmin           = nmin,
    asSF           = TRUE)

  message(blue("\n\n*****  ASSIGNING GRID IDENTIFIERS *****\n"))
  
  spDataGridDF <- prepSpDataWithGridID(
    sfSpeciesDF    = sfSpeciesDF,
    sfGrid         = sfGrid,
    removeDups     = removeDups,
    filterByMinObs = filterByMinObs,
    nmin           = nmin,
    getCounts      = FALSE)
                          
  message(blue("\n\n***** READING OCCURRENCES' ENVIRONMENTAL DATA PER YEAR *****\n\n"))
  
  envDataOcc <- getEnvDataByYear(
    dfPathList  = dfPathList,
    spDataGrid  = spDataGridDF,
    progressBar = progressBar)
  
  
  if(!is.null(spName)){
    
    message(blue("\n\n***** GENERATING PSEUDO-ABSENCES FOR TARGET SPECIES *****\n\n"))
    
    PseudoAbsData <- doRandomPA(
      spName      = spName,
      dfPathList  = dfPathList,
      envDataOcc  = envDataOcc,
      nPAsets     = nPAsets,
      nPAperSet   = nPAperSet,
      progressBar = progressBar)
    
    message(blue("\n\n***** EXPORTING DATA FOR TARGET SPECIES *****\n\n"))
    
    return(list(trainData = bind_rows(envDataOcc %>% filter(SpeciesName == spName),
                                      mergeDFsInList(PseudoAbsData)),
                trainMatrix = doTrainMatrix(PseudoAbsData)
    ))

  }else{
    
    message(blue("\n\n***** GENERATING PSEUDO-ABSENCES BY SPECIES *****\n\n"))

    spNames <- getSpeciesNames(spDataGridDF)
    
    if(progressBar){
      pb <- txtProgressBar(1, length(spNames), style=3)
    } 
    
    i<-0
    
    for(spName in spNames){
    
      # Increment counter  
      i<-i+1

      message(green("\nGenerating pseudo-absences data for",spName,"by year:\n"))
      
      # Generate the pseudo-absences per species
      PseudoAbsData <- doRandomPA(
        spName      = spName,
        dfPathList  = dfPathList,
        envDataOcc  = envDataOcc,
        nPAsets     = nPAsets,
        nPAperSet   = nPAperSet,
        progressBar = progressBar)
      
      trainData <- bind_rows(envDataOcc %>% filter(SpeciesName == spName),
                            mergeDFsInList(PseudoAbsData))
      
      trainMatrix <- doTrainMatrix(PseudoAbsData)
      
      write_rds(trainData, 
                paste(outDir,"/",gsub("\\ +","_",spName),
                      "_TrainData_PPA.rds",sep=""))
      write_rds(trainMatrix, 
                paste(outDir,"/",gsub("\\ +","_",spName),
                      "_TrainMatrix_PPA.rds",sep=""))
      write_csv(trainData, 
                paste(outDir,"/",gsub("\\ +","_",spName),
                      "_TrainData_PPA.csv",sep=""))
      write_csv(as.data.frame(trainMatrix),
                paste(outDir,"/",gsub("\\ +","_",spName),
                      "_TrainMatrix_PPA.csv",sep=""))
     
      message(yellow("\n\nFinished species:",spName,"\n\n"))
      setTxtProgressBar(pb, i)                            
    }
    return(envDataOcc)
  }
}



#' Perform iterative variable elimination and selection
#' 
#' This function combines Random Forests variable importance estimation and pairwise correlation 
#' analysis to first rank variables according to their predictive importance and then it performs 
#' correlation to iteratively select the less correlated ones. Variable selection procedures can 
#' be applied independently by groups of predictors. It also performs preliminary steps based on 
#' caret package near-zero variance function and correlatio elimination using a very-high threshold 
#' value.
#' 
#' 
#' @param trainData A train data dataframe object created by \code{createTrainData}
#' @param trainMatrix A train matrix object created by \code{createTrainData} displaying 
#' which rows to use in each set of pseudo-absences round
#' @param predVars Character vector with the names of predictor variables
#' @param ntree Number of trees to grow in Random Forests algorithm (default: 1000)
#' @param nmax Maximum number of predictors to select. If set to "10-rule" it uses Harrell's thumbrule of 
#' 10 observations per predictor (default: "10-rule")
#' @param tenRuleNmax If applying Harrell's 10-rule this parameter defines a maximum number of predictors 
#' to avoid using a very large of factors (default: NULL, i.e. not used)
#' @param thresh Correlation threshold used in iterative variable selection. If above this value the 
#' predictor will get discarded (default: |thresh| = 0.8)
#' @param aggFun Aggregation function used to merge values from all PA sets (default: median)
#' @param method Correlation method used (default: "spearman")
#' @param varGroups Used to perform the variable selection by group with an equal number of predictors 
#' for each group. The input should be a table with at least two columns identifying the name of the 
#' predictor variable (column: "VarName") and other one with the group name (column: "VarGroup"). 
#' Column names must be used as these ones
#' @param doNZV Perform caret package near-zero variance preliminary analysis? (default: TRUE)
#' @param doHighCor Perform caret package preliminary correlation analyses based on a very-high correlation
#' value? (default: TRUE)
#' @param highCorVal Value used to perfrom the preliminary caret correlation elimination 
#' (default: |highCorVal| = 0.95). 
#' @param pbar Print progress bar? (default: TRUE)
#' @param verbose Print messages? (default: TRUE)
#' 
#' @return A character vector with selected variables
#' 
#' @note The near zero variance analyses and the "very-high correlation" filtering/removal will be 
#' perform previously to the RF models and the iterative selection processes. This can be considered as 
#' a first screening.
#' 

multiRoundIterVarSel <- function(trainData, 
                                 trainMatrix, 
                                 predVars, 
                                 ntree        = 1000,
                                 nmax         = "10-rule",
                                 tenRuleNmax  = NULL,  
                                 thresh       = 0.8, 
                                 aggFun       = median,
                                 method       = "spearman",
                                 varGroups    = NULL, 
                                 doNZV        = TRUE, 
                                 doHighCor    = TRUE, 
                                 highCorVal   = 0.95,
                                 pbar         = TRUE, 
                                 verbose      = TRUE,
                                  ...){
  
  if(verbose) message(blue("\n\n***** DOING PRELIMINARY STUFF *****\n\n"))
  
  x <- trainData[,predVars]
  
  if(doNZV){
    nzvRem <- caret::nearZeroVar(x[,predVars])
    x <- select(x, -all_of(nzvRem))
  }
  if(doHighCor){
    cmat <- cor(x, method=method)
    corRem <- caret::findCorrelation(cmat, cutoff = highCorVal)
    x <- select(x, -all_of(corRem))
  }

  
  if(verbose) message(blue("\n\n***** RUNNING RANDOM FOREST BY PA SET *****\n\n"))
  
  
  if(pbar) pb <- txtProgressBar(min = 1, max = ncol(trainMatrix), style = 3)
  
  
  for(i in 1:ncol(trainMatrix)){

    trainVec <- trainMatrix[,i]
    y.train <- as.factor(trainData[trainVec,"pa"])
    x.train <- x[trainVec, ]
    
    if(nmax == "10-rule"){
      
      nmax <- round(nrow(x.train) / 10)
      
      if(!is.null(tenRuleNmax)){
        if(nmax > tenRuleNmax){
          nmax <- tenRuleNmax
        }
      }
      
      if(!is.null(varGroups)){
        nGroups <- length(unique(varGroups$VarGroup))
        nmax <- round(nmax / nGroups)
        groupNames <- unique(unlist(varGroups[,"VarGroup"]))
      }else{
        nGroups <- 1
      }
    }
    
    rf <- randomForest(y = y.train, x = x.train, ntree=ntree)
    
    tmpImp <- importance(rf)
    tmpImp <- data.frame(VarName = rownames(tmpImp), 
                         importance = tmpImp[,1], stringsAsFactors = FALSE) %>% 
      arrange(VarName)
    colnames(tmpImp) <- c("VarName",paste("r",i,sep="_"))
    
    if(i==1){
      impMat <- tmpImp
    }else{
      impMat <- bind_cols(impMat,
                          tmpImp %>% dplyr::select(2))
    }
    
    setTxtProgressBar(pb, i)
    
  }
  
  
  if(verbose) message(blue("\n\n***** ITERATIVE SELECTION USING IMPORTANCE & CORRELATION *****\n\n"))
  
  
  vimpDF <- data.frame( VarName = impMat[,"VarName"],
                        vimpAvg  = apply(impMat[,-1], 1, FUN = aggFun),
                        vimpStd  = apply(impMat[,-1], 1, FUN = sd),
                        vimpMAD  = apply(impMat[,-1], 1, FUN = mad)) %>% 
    arrange(desc(vimpAvg))
  
  if(!is.null(varGroups)){
    
    vimpDF <- vimpDF %>% left_join(varGroups, by ="VarName")
    vimpDF_init <- vimpDF
  }
  
  

  for(j in 1:nGroups){
    
    if(!is.null(varGroups)){
      vimpDF <- vimpDF_init %>% filter(VarGroup == groupNames[j])
      message(yellow("\nSelecting in group:",groupNames[j],"\n"))
    }
    
    selVars <- vimpDF[1,1]
    
    for(i in 2:nrow(vimpDF)){
      
      activeVar <- vimpDF[i,1]
      cmat <- cor(trainData[,c(selVars,activeVar)], method = method)
      corVec <- abs(cmat[lower.tri(cmat)])
      
      if(sum(corVec > thresh) == 0){
        selVars <- c(selVars,activeVar)
      }
      
      if(!is.null(nmax)){
        
        if(length(selVars) >= nmax){
          if(verbose) message(yellow("\nReached the maximum number of variables (nmax)!\n"))
          break
        }
      }
    }
    
    if(j==1){
      selVarsAll <- c(selVars)
    }else{
      selVarsAll <- c(selVarsAll,selVars)
    }
  }
  
  if(!is.null(varGroups)){
    attr(selVarsAll,"vimpDF") <- vimpDF_init
  }else{
    attr(selVarsAll,"vimpDF") <- vimpDF
  }
  
  return(selVarsAll)
}


## Functions for post-processing modelling HS time series ---------------------------------


#' Majority filter
#' 
#' Check if a value is in majority
#' 
#' @param x Value count
#' @param n Total of cases/observations
#' 
#' @return A boolean value with 1 if is in majority and 0 otherwise

majFilter <- function(x, n){
  if(is.nan(x) || is.na(x)){
    return(NA)
  }else{
    if(x >= n/2){
      return(1)
    }else{
      return(0)
    }
  }
} 


#' @rdname majFilter

MajFilter <- function(x, n){ 
  sapply(x, majFilter, n=n)  
}


#' Simple plot of habitat change parameters
#' 
#' Plots the amount of adequate habitat through time
#' 
#' @param y Data to plot in y
#' @param x Years
#' @param ... Additional params passed to plot
#' 
#' @return A plot

plotHabChangePars <- function(y, x= seq(1985, 2019, by=2),...){
  plot(x=x, y=y, ...)
  lines(x=x, y=y)
  abline(v=x, lty="dashed", col="light grey")
}



## Functions to calculate the Theil-Sen trend slope and the p-value ----------------------


#' Theil-Sen slope magnitude value
#' 
#' Calculates Theil-Sen slope magnitude value
#' 
#' @param x A numeric vector with a regular time series
#' @param ... Other parameters passed to \code{sens.slope} function
#' 
#' @return A numeric value with the slope magitude

senSlope <- function(x,...) sens.slope(x,...)$estimates


#' Theil-Sen trend p-value 
#' 
#' Calculates Theil-Sen trend slope p-value
#' 
#' @param x A numeric vector with a regular time series
#' @param ... Other parameters passed to \code{sens.slope} function
#' 
#' @return A numeric value with the slope magitude

senSlopePval <- function(x,...) sens.slope(x,...)$p.value


#' Theil-Sen trend 
#' 
#' Calculates Theil-Sen slope magnitude value
#' 
#' @param x A numeric vector with a regular time series
#' @param ... Other parameters passed to \code{sens.slope} function
#' 
#' @return A numeric value with the slope magitude

senTrendSlope <- function(x, na.rm = TRUE, ...){
  
  if(all(is.na(x) | is.nan(x))){
    return(c(NA,NA))
  }else{
    if(na.rm){
      x <- na.omit(x)
    }
  }
  x <- as.numeric(x)
  sens <- try(suppressWarnings(sens.slope(x,...)))
  
  if(inherits(sens,"try-error")){
    return(c(NA,NA))
  }else{
    return(as.numeric(c(sens[["estimates"]], sens[["p.value"]])))
  }
}


