


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


grid1k_wgs84 <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid_WGS84.shp")

dfPathList <- list.files("./DATA_/TABLES/_PRED_DATASETS/", ".rds$", full.names = TRUE)

spData <- readxl::read_excel("./DATA_/TABLES/SpeciesData/Amphibians/endemic_amphibia_Espinhaco_v2-1.xlsx")

#length(unique(spData$Species))


## -------------------------------------------------------------- ##

getYears <- function(x) return(x %>% pull(Year) %>% unique %>% sort)

getSpeciesNames <- function(x) return(x %>% pull(SpeciesName) %>% unique %>% sort)

getIDs <- function(x) return(x %>% pull(ID))

getYearIDs <- function(x, yr) return(x %>% filter(Year == yr) %>% pull(ID))

readDatasetByYear <- function(fileList, year){
  fpath <- fileList[grepl(year,fileList)]
  return(read_rds(fpath))
}

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
      filter (!duplicated(cd))
      #distinct(SpeciesName, ID, Year, 
      #                         .keep_all = TRUE)
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



getEnvDataByYear <- function(dfPathList,
                             spDataGrid,
                             progressBar = TRUE){

  if("spData" %in% names(spDataGrid)){
    spDataGrid <- spDataGrid[["spData"]]
  }
  
  if(inherits(spDataGrid,"sf")){
    spDataGrid <- spDataGrid %>% st_drop_geometry()
  }
  
  #print(spDataGrid)
  
  spDataGrid <- data.frame(pa = 1, as.data.frame(spDataGrid))
  
  yrs <- getYears(spDataGrid)
  
  if(progressBar){
    pb <- txtProgressBar(1,length(yrs),style = 3)
    cat(green("\nReading presence data by year:\n"))
  }
  
  for(i in 1:length(yrs)){
    
    spDataGridTmp <- spDataGrid %>% filter(Year == yrs[i])
    envDataYr <- readDatasetByYear(dfPathList, yrs[i]) #%>% 
                  #filter(ID %in% getIDs(spDataGridTmp))
    
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
    cat(green("\nGenerating pseudo-absences data by year:\n"))
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
  
  
  cat(blue("\n\n--- PREPARING AND FILTERING DATA ---\n\n"))
  
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

  cat(blue("\n\n--- ASSIGNING GRID IDENTIFIERS ---\n\n"))
  
  spDataGridDF <- prepSpDataWithGridID(
    sfSpeciesDF    = sfSpeciesDF,
    sfGrid         = sfGrid,
    removeDups     = removeDups,
    filterByMinObs = filterByMinObs,
    nmin           = nmin,
    getCounts      = FALSE)
  
  # spDataGridDF %>% 
  #   mutate(cd=paste(gsub("\\ +","_",SpeciesName),ID,Year,sep="_")) %>% 
  #   filter (!duplicated(cd)) %>% 
  #   #distinct(cd, .keep_all = TRUE) %>% 
  #   arrange(SpeciesName,ID) %>% 
  #   View
                          
  cat(blue("\n\n--- READING OCCURRENCES' ENVIRONMENTAL DATA PER YEAR ---\n\n"))
  
  envDataOcc <- getEnvDataByYear(
    dfPathList  = dfPathList,
    spDataGrid  = spDataGridDF,
    progressBar = progressBar)
  
  
  if(!is.null(spName)){
    
    cat(blue("\n\n--- GENERATING PSEUDO-ABSENCES FOR TARGET SPECIES ---\n\n"))
    
    PseudoAbsData <- doRandomPA(
      spName      = spName,
      dfPathList  = dfPathList,
      envDataOcc  = envDataOcc,
      nPAsets     = nPAsets,
      nPAperSet   = nPAperSet,
      progressBar = progressBar)
    
    cat(blue("\n\n--- EXPORTING DATA FOR TARGET SPECIES---\n\n"))
    
    return(list(trainData = bind_rows(envDataOcc %>% filter(SpeciesName == spName),
                                      mergeDFsInList(PseudoAbsData)),
                trainMatrix = doTrainMatrix(PseudoAbsData)
    ))

  }else{
    
    cat(blue("\n\n--- GENERATING PSEUDO-ABSENCES BY SPECIES ---\n\n"))

    spNames <- getSpeciesNames(spDataGridDF)
    
    if(progressBar){
      #cat(green("\nGenerating pseudo-absences data by species/year:\n"))
      pb <- txtProgressBar(1, length(spNames), style=3)
    } 
    
    i<-0
    for(spName in spNames){
      i<-i+1

      PseudoAbsData <- doRandomPA(
        spName      = spName,
        dfPathList  = dfPathList,
        envDataOcc  = envDataOcc,
        nPAsets     = nPAsets,
        nPAperSet   = nPAperSet,
        progressBar = progressBar)
      
      #cat(blue("\n\n--- EXPORTING DATA FOR SPECIES:",spName," ---\n\n"))
      
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
     
      cat(yellow("\n\nFinished species:",spName,"\n\n"))
      setTxtProgressBar(pb, i)                            
    }
    return(envDataOcc)
  }
}


## -------------------------------------------------------------- ##





spDatasf <- prepSpData(spData, filterByMinObs = FALSE)

spDataGrid <- prepSpDataWithGridID(spDatasf, grid1k_wgs84, 
                         filterByMinObs = FALSE, getCounts = TRUE)

View(spDataGrid[["spCounts"]])
#View(spDataGrid[["spCounts"]])


envDataOcc <- getEnvDataByYear(dfPathList,
                             spDataGrid)


trainDataFull <- 
createTrainData (spData = spData, 
                 spName = "Boana botumirim",
                 sfGrid = grid1k_wgs84,
                 dfPathList = dfPathList,
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
                 progressBar    = TRUE
)



spData = spData 
spName = "Boana botumirim"
sfGrid = grid1k_wgs84
dfPathList = dfPathList
spNameCol      = "Species" 
lonCol         = "Lon" 
latCol         = "Lat" 
yearCol        = "Year"
yrStart        = 1985 
yrEnd          = 2019 
filterByMinObs = TRUE
nmin           = 30
removeDups     = TRUE
nPAsets        = 10
nPAperSet      = "equal"
progressBar    = TRUE


createTrainData (spData           = spData, 
                   spName         = NULL,
                   sfGrid         = grid1k_wgs84,
                   dfPathList     = dfPathList,
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
                   outDir         = "./DATA_/TABLES/_TRAIN_DATASETS/")




