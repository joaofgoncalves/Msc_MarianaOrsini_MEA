


library(terra)
library(sf)
library(dplyr)
library(ggplot2)

rm(list=ls())

# plotOutDir <- "./OUT/MODS/R02_Plots"
#plotOutDir <- "./OUT/MODS/R02_Plots_CropSE"
plotOutDir <- "./OUT/MODS/R02_Plots_CropRBSE"

# spGroups <- list.dirs("./OUT/MODS/R02", recursive = FALSE)
#spGroups <- list.dirs("./OUT/MODS/R02_CropSE", recursive = FALSE)
spGroups <- list.dirs("./OUT/MODS/R02_CropRBSE", recursive = FALSE)


## ---------------------------------------------------------------------------------- ##

yrs_ <- seq(1985,2019,by=2)
yrs <- rep(seq(1985,2019,by=2),each=2)

selVarsDF <- readxl::read_excel("./DATA_/TABLES/LULC_Habitats_Species-simple-v2.xlsx")



k <- 0
for(spGroup in spGroups){
  
  spNames <- list.dirs(spGroup, recursive = FALSE) 
  
  spGroup_ <- basename(spGroup)
  
  
  for(spName in spNames){
    
    k <- k + 1
    spName_ <- basename(spName)
    
    spNameFull <- selVarsDF$Species[selVarsDF$SpeciesComp %in% spName_]
    
    
    spFileList <- list.files(spName, pattern=".tif$", full.names = TRUE)
    
    binFiles <- spFileList[grepl("TSSbin",spFileList)] 
    #hsbFiles <- spFileList[!grepl("TSSbin",spFileList)]
    
    
    rstBin <- rast(binFiles)
    
    freqTableSpInit <- terra::freq(rstBin) %>% 
      as.data.frame() 
    
    
    if(!all(freqTableSpInit$value == 0)){
      
      areaTotal <- sum(freqTableSpInit[1:2,"count"])
      
      freqTableSpInit <- freqTableSpInit %>% 
        mutate(yrs = yrs, spName = spName_, spGroup = spGroup_) %>% 
        mutate(percArea = (count/areaTotal) * 100) %>% 
        select(spGroup, spName, yrs,value,count, percArea)
      
      freqTableSp <- freqTableSpInit %>% 
        filter(value == 1) %>% 
        mutate(yrs = yrs_, spName = spName_, spGroup = spGroup_) %>% 
        select(spGroup, spName, yrs,value,count, percArea)
      
      
      if(k==1 | !exists("freqTableAllSp")){
        freqTableAllSp <- freqTableSpInit
      }else{
        freqTableAllSp <- bind_rows(freqTableAllSp, freqTableSpInit)
      }
      
      
    }else{
      
      freqTableSpInit <- freqTableSpInit %>% 
        mutate(yrs = yrs_, spName = spName_, spGroup = spGroup_) %>% 
        mutate(percArea = 100) %>% 
        select(spGroup, spName, yrs,value,count, percArea)
      
      
      if(k==1 | !exists("freqTableAllSp")){
        freqTableAllSp <- freqTableSpInit
      }else{
        freqTableAllSp <- bind_rows(freqTableAllSp, freqTableSpInit)
      }
      
      
      cat(spNameFull,"has no suitable locations!!\n\n")
      next
    }
    
    
    mn <- mean(freqTableSp$count)
    
    g <- ggplot(freqTableSp,aes(x=yrs,y=count))+ 
      geom_line() + 
      geom_smooth(se = FALSE) + 
      geom_point() + 
      labs(title = paste(spGroup_,"|",spNameFull)) +
      xlab("Year") + 
      ylab(bquote("Suitable habitat ("~Km^2~")")) +
      theme_bw()
    
    g1 <- ggplot(freqTableSp,aes(x=yrs,y=((count - mn) / mn)*100))+ 
      geom_line() + 
      geom_smooth(se = FALSE) + 
      geom_point() + 
      labs(title = paste(spGroup_,"|",spNameFull)) +
      xlab("Year") + 
      ylab(bquote("Suitable habitat (% of avg.)")) +
      theme_bw()
    
    #plot(g)
    
    fn <- paste(plotOutDir,"/SuitableAreaTrend_",spGroup_,"_",spName_,"_v1.png",sep="")
    
    ggsave(filename = fn, plot = g, width = 6, height = 5)
  
    
    fn <- paste(plotOutDir,"/SuitableAreaTrend_",spGroup_,"_",spName_,"_relChange_v1.png",sep="")
    
    ggsave(filename = fn, plot = g1, width = 6, height = 5)
    
    cat("Finished...",spNameFull,"\n\n\n")  
  }
  
}

#write.csv(freqTableAllSp, paste(plotOutDir,"/_freqTableAllSpByYear_CropSE.csv",sep=""), row.names = FALSE)


write.csv(freqTableAllSp, paste(plotOutDir,"/_freqTableAllSpByYear_CropRBSE.csv",sep=""), row.names = FALSE)


