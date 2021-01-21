

library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(fasterize)
library(exactextractr)

st_erase <- function(x, y) st_difference(x, st_union(st_combine(y)))

calcChangeFraction <- function(x, cf = NULL){
  
  OUT <- matrix(NA, nrow = length(x), ncol=2)
  colnames(OUT) <- c("PID","PropChange")
  
  for(i in 1:length(x)){
    
    DF <- x[[i]]

    if(is.null(cf)){
      TMP <- DF %>% 
        group_by(value) %>% 
        summarize(nPix = sum(coverage_fraction), .groups ="drop_last")
    }else{
      TMP <- DF %>% 
        filter(coverage_fraction >= cf) %>% 
        group_by(value) %>% 
        summarize(nPix = n(), .groups ="drop_last")
    }
    
    if(nrow(TMP)==1){
      if(TMP[1,1] == 0){
        OUT[i,] <- c(i, 0)
      }else if(TMP[1,1] == 1){
        OUT[i,] <- c(i, 1)
      }
    }else{ 
      OUT[i,] <- c(i, as.numeric(TMP[2,2] / sum(TMP[,2])))
    }
  }
  return(as.data.frame(OUT))
}

## ----------------------------------------------------------------------------- ##

r_2000 <- c("./DATA_/RASTER/LULC_MG_UTM_23S/mapbiomas-brazil-collection-50-minasgerais-2000_WGS84_UTM23S.tif")
r_2019 <- c("./DATA_/RASTER/LULC_MG_UTM_23S/mapbiomas-brazil-collection-50-minasgerais-2019_WGS84_UTM23S.tif")

r2000 <- raster(r_2000)
r2019 <- raster(r_2019)

#rDiff <- r2000 != r2019
#writeRaster(rDiff,"./DATA_/RASTER/LULC_MG_UTM_23S/MapBiomas_MG_Diff_2000_2019.tif")
rDiff <- raster("./DATA_/RASTER/LULC_MG_UTM_23S/MapBiomas_MG_Diff_2000_2019.tif")

shpMinasGerais <- read_sf("./DATA_/VECTOR/gadm36_BRA_shp/gadm36_BRA_1.shp") %>% 
  filter(NAME_1 == "Minas Gerais") %>% 
  st_transform("EPSG:32723") %>% 
  mutate(val=0)

shpPath <- "./DATA_/VECTOR/MG_minas_jazidas/MG_distribuicao_de_jazidas_e_minas.shp"

minasJazidas <- read_sf(shpPath)

minasJazidas23S <- st_transform(minasJazidas, "EPSG:32723")


## ----------------------------------------------------------------------------- ##


dists <- c(150, 300, 600, 1000, 2000, 5000)
N <- 50

outDiffs <- list()
outPvals <- list()


for(buffDist in dists){
  
  
  cat("\n\nRunning buffer distance =", buffDist, ".......\n\n")
  
  minasJazidas23S_buff <- minasJazidas23S %>% 
    st_buffer(dist = buffDist) %>% 
    st_union() %>% 
    st_cast(to="POLYGON") %>% 
    st_as_sf()
  
  extr_Diff_Minas <- exact_extract(rDiff, minasJazidas23S_buff, fun = NULL, 
                                 force_df=TRUE, progress=FALSE)
  
  changeFractions_Minas <- calcChangeFraction(extr_Diff_Minas)

  avgChangeDiffs <- vector(mode="numeric", length=N)
  pvals <- vector(mode="numeric", length=N)
  
  pb <- txtProgressBar(max=N, style=3)  
  
  for(i in 1:N){
    
    randPts_buff <- st_sample(shpMinasGerais, 
                              size = nrow(minasJazidas23S), 
                              type = "random") %>% 
      st_buffer(dist = buffDist) %>% 
      st_union() %>% 
      st_cast(to="POLYGON") %>% 
      st_as_sf() %>% 
      st_erase(minasJazidas23S_buff) %>% 
      st_cast(to="POLYGON")
    
    extr_Diff_rnd <- exact_extract(rDiff, randPts_buff, fun = NULL, 
                                   force_df=TRUE, progress=FALSE)
    
    changeFractions_rnd <- calcChangeFraction(extr_Diff_rnd)
    
    avgChangeDiffs[i] <- mean(changeFractions_Minas[,2]) - mean(changeFractions_rnd[,2])
    
    wt <- wilcox.test(x = changeFractions_Minas[,2],
                      y = changeFractions_rnd[,2],
                      alternative = "greater")
    
    pvals[i] <- wt$p.value
    
    setTxtProgressBar(pb,i)
  }
  
  outDiffs[[paste("d=",buffDist,sep="")]] <- avgChangeDiffs
  outPvals[[paste("d=",buffDist,sep="")]] <- pvals
  
}


## ----------------------------------------------------------------------------- ##


for(i in 1:length(outDiffs)){
  
  TMP <- data.frame(dist = dists[i],
                    avgChanDiffs = outDiffs[[i]],
                    pVals = outPvals[[i]],
                    pVals_adjBY = p.adjust(outPvals[[i]], method = "BY"))
  
  if(i==1){
    outDiffsDF <- TMP
  }else{
    outDiffsDF <- bind_rows(outDiffsDF, TMP)
  }
}

outDiffsDF %>% 
  group_by(dist) %>% 
  summarize_all(.funs = list(avg = mean))


