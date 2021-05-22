

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(readr)

#fl <- list.files("./DATA_/TABLES/VARS_PREP/LULC/", pattern=".csv$",full.names = TRUE)

# fl <- list.files("E:/Projects/MSc_MarianaOrsini/DATA/TABLES/VARS_PREP/EVI2", 
#                  pattern=".csv$",full.names = TRUE)

fl <- list.files("../../../GEE/", 
                 pattern=".csv$",full.names = TRUE)


for(i in 1:length(fl)){
  
  yr <- rep((1985:2019),2)[i]
  
  #tmpDF <- read.csv(fl[i])
  tmpDF <- read_csv(fl[i])
  

  
  if(grepl("_avg_",fl[i])){
    
    tmpDF <- tmpDF %>% 
      select(ID, `mean`) %>% 
      mutate(ID = as.integer(ID)) %>% 
      mutate(`mean` = round(`mean` * 10000)) %>% 
      arrange(ID)
    
    tmpDF <- tmpDF %>%
      #`colnames<-`(c("ID",paste("EVI2med_avg_",yr,sep="")))
      `colnames<-`(c("ID",paste("NDWImed_avg_",yr,sep="")))
  }else{
    
    tmpDF <- tmpDF %>% 
      select(ID, `stdDev`) %>% 
      mutate(ID = as.integer(ID)) %>% 
      mutate(`stdDev` = round(`stdDev` * 10000)) %>% 
      arrange(ID)
    
    tmpDF <- tmpDF %>%
    #`colnames<-`(c("ID",paste("EVI2med_std_",yr,sep=""))) 
    `colnames<-`(c("ID",paste("NDWImed_std_",yr,sep=""))) 
  }

  
  if(i==1){
    EVI2DF <- tmpDF
  }else{
    EVI2DF <- EVI2DF %>% left_join(tmpDF, by = "ID")
  }
  
  cat("\n\n-> Processing year:",yr,"\n\n")
}

#write_csv(EVI2DF,"./DATA_/TABLES/VARS/EVI2DF_ByYear1985_2019-v1.csv")
#write_rds(EVI2DF,"./DATA_/TABLES/VARS/EVI2DF_ByYear1985_2019-v1.rds")

write_csv(EVI2DF,"./DATA_/TABLES/VARS/NDWIDF_ByYear1985_2019-v1.csv")
write_rds(EVI2DF,"./DATA_/TABLES/VARS/NDWIDF_ByYear1985_2019-v1.rds")


rr <- sample(1:nrow(EVI2DF),1)
par(mfcol=c(2,1))
plot(x= 1985:2019, y=EVI2DF[rr,2:36])
lines(x= 1985:2019, y=EVI2DF[rr,2:36])
plot(x= 1985:2019, y=EVI2DF[rr,37:71])
lines(x= 1985:2019, y=EVI2DF[rr,37:71])



## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##


fl <- list.files("./DATA_/TABLES/VARS_PREP/LULC/", pattern=".csv$",full.names = TRUE)

newClasses <- readxl::read_excel("./DATA_/TABLES/MapBiomasNewClasses-v1.xlsx")


for(i in 1:length(fl)){
  
  yr <- (1985:2019)[i]
  
  #tmpDF <- read.csv(fl[i])
  tmpDF <- read_csv(fl[i])
  
  idx <- as.integer(gsub("cl_","",colnames(tmpDF)[-1]))
  clCodesToUse <- newClasses$ClassCode[idx]
  
  tmpDFrelFreq <- tmpDF[,-1] / apply(tmpDF[,-1],1,sum)
  colnames(tmpDFrelFreq) <- paste(clCodesToUse,yr,sep="_")
  
  tmpDF <- bind_cols(tmpDF[,1], tmpDFrelFreq)
  
  if(i==1){
    LULCDF <- tmpDF
  }else{
    LULCDF <- LULCDF %>% left_join(tmpDF, by = "ID")
  }
  
  cat("\n\n-> Processing year:",yr,"\n\n")
}

write_csv(LULCDF,"./DATA_/TABLES/VARS/LULCDF_ByYear1985_2019-v1.csv")
write_rds(LULCDF,"./DATA_/TABLES/VARS/LULCDF_ByYear1985_2019-v1.rds")

# LULCDF %>% select(ends_with("_1985"))

