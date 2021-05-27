
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
library(trend)
library(diptest)

## -------------------------------------------------------------------------------------- ##

#setwd("C:/MyFiles/R-dev/Msc_MarianaOrsini_MEA")

# Read ancillary functions
source("./RCODE/_LibMainFunctions-v1.R")

## -------------------------------------------------------------------------------------- ##


# Read 1-Km grid data in vector and raster
grid1k <- read_sf("./DATA_/VECTOR/Area_estudo/espinhaco_grid.shp")
rstGrid1K <- raster("./DATA_/RASTER/Grid1K_IDs_WGS84UTM23S_v1.tif")

nYears <- 18



## -------------------------------------------------------------------------------------- ##

#i=1
#dirList <- list.dirs("./OUT/MODS/R01", recursive = FALSE)
dirList <- list.dirs("C:/Users/JG/Google Drive (up200103770@g.uporto.pt)/CIBIO_Desktop/R_projects/TeseMariana_CIBIOPC/OUT/MODS_CIBIOPC/R01", 
                     recursive = FALSE)


for(i in 1:length(dirList)){
  
  # Read data from the main rds file for the species
  # The rds file contains HS and binary projections for all dates 
  # and all years in the time series
  # Binary projections column name suffix is: _TSSbin
  # 
  fpath <- list.files(dirList[i], full.names = TRUE, pattern=".rds$")
  fpath <- fpath[grepl("_allProjections",fpath)]
  allProjs <- read_rds(fpath)
  
  # Calculate the row-wise sum of predictions across all years
  # Each row represents a 1 Km-square quadrat
  #
  predSum <- allProjs %>% 
    select(ends_with("_TSSbin")) %>% 
    apply(MARGIN = 1, FUN = sum) %>% 
    data.frame(PredSum =.) %>% 
    mutate(PredMaj = MajFilter(PredSum, nYears)) %>% 
    bind_cols(allProjs %>% select(ID)) %>% 
    select(ID,PredSum,PredMaj)
  
  # Calculate the annual totals of suitable locations (i.e., available suitable habitat)
  #
  predsByYear <- allProjs %>% 
    select(ends_with("_TSSbin")) %>% 
    apply(MARGIN = 2, FUN = sum, na.rm=TRUE)
  
  # Calculate median habitat suitability for all suitable 
  # pixels (across time) found in the time series
  #
  medHsByYear <- allProjs %>% 
    select(ends_with("_HS_ens")) %>% 
    filter(predSum$PredSum > 0) %>% 
    apply(MARGIN = 2, FUN = median, na.rm=TRUE)
  
  plotHabChangePars(y=predsByYear, xlab="Year",ylab="Suitable Area (km2)", 
                    main=basename(dirList[i]))
    
  plotHabChangePars(y=medHsByYear, xlab="Year", 
                    ylab="Median habitat suitability (all suitable locations)", 
                    main=basename(dirList[i]))
  
  
  # Calculate Theil-Sen slope monotonic trend -------------------------------------
  
  
  # Calculate Sen's slope per row/ 1Km-quadrat
  senSlopeTrend <- allProjs %>% 
    select(ends_with("_HS_ens")) %>% 
    apply(MARGIN = 1, FUN = senTrendSlope) %>% 
    t %>% as.data.frame %>% 
    `colnames<-`(c("SenSlope","Pvalue"))
  
  senSlopeTrend <- bind_cols(allProjs %>% select(ID),
                      senSlopeTrend)
  # Merge data
  grid1KPredSumTrends <- grid1k %>% 
    left_join(predSum, by="ID") %>% 
    left_join(senSlopeTrend, by="ID")
  
  
  ## Prepare masks ---------------------------------------------------------------
  
  # Mask with pixels that were suitable the majority of years
  maskMaj <- fasterize::fasterize(grid1KPredSumTrends, raster = rstGrid1K, field = "PredMaj")
  # Mask with the number of time a pixels is suitable
  rstPredSum <- fasterize::fasterize(grid1KPredSumTrends, raster = rstGrid1K, field = "PredSum")
  # Mask with pixels that are suitable in any year of the series
  maskAny <- rstPredSum > 0
  
  writeRaster(maskMaj, paste(dirList[i],"/",basename(dirList[i]),"_PredMajority.tif",sep=""))
  writeRaster(rstPredSum, paste(dirList[i],"/",basename(dirList[i]),"_PredSum.tif",sep=""))
  
  ## Trend Sen-Theil slope maps -------------------------
  
  tmpRst <- fasterize::fasterize(grid1KPredSumTrends, raster = rstGrid1K, field = "SenSlope")
  writeRaster(tmpRst, paste(dirList[i],"/",basename(dirList[i]),"_TrendSenSlope_AllPixs.tif",sep=""))
  
  tmpRstMaj <- tmpRst
  tmpRstMaj[maskMaj == 0] <- NA
  writeRaster(tmpRstMaj, paste(dirList[i],"/",
                            basename(dirList[i]),"_TrendSenSlope_MajPixs.tif",sep=""))
  
  tmpRstAny <- tmpRst
  tmpRstAny[maskAny == 0] <- NA
  writeRaster(tmpRstAny, paste(dirList[i],"/",
                            basename(dirList[i]),"_TrendSenSlope_AnyPixs.tif",sep=""))
  
  ## P-value maps -----------------------------------------------------------------
  
  tmpRst <- fasterize::fasterize(grid1KPredSumTrends, raster = rstGrid1K, field = "Pvalue")
  writeRaster(tmpRst,paste(dirList[i],"/",basename(dirList[i]),"_TrendSenPvalue_AllPixs.tif",sep=""))
  
  tmpRstMaj <- tmpRst
  tmpRstMaj[maskMaj == 0] <- NA
  writeRaster(tmpRstMaj, paste(dirList[i],"/",
                            basename(dirList[i]),"_TrendPvalue_MajPixs.tif",sep=""))
  
  tmpRstAny <- tmpRst
  tmpRstAny[maskAny == 0] <- NA
  writeRaster(tmpRstAny, paste(dirList[i],"/",
                            basename(dirList[i]),"_TrendPvalue_AnyPixs.tif",sep=""))
  
}

