

library(terra)
library(rasterVis)
library(raster)

# The aggregation factor
aggFactor <- 167


# Read the transition map
# DON'T FORGET TO CHANGE THE INPUT FILE
rst <- rast("C:/Users/JG/Desktop/rbse1_transitionsmap_1985_2019_v1.tif")


# The value cl in the loop equals each transition passed to the binary raster test

for(cl in 1:13){
  
  # Create a raster binary layer for each transition class
  rstBin <- (rst == cl)
  
  # Perform raster aggregation by a scale factor defined in aggFactor
  aggRst_Cl <- terra::aggregate(rstBin, fact= aggFactor, fun="sum")
  
  # Calculate the % of each transition class
  aggRst_Cl <- (aggRst_Cl / aggFactor^2) * 100
  
  # Write the data to file
  # DON'T FORGET TO CHANGE THE OUPUT DIRECTORY
  terra::writeRaster(aggRst_Cl, paste("C:/Users/JG/Desktop/maps/TransitionPerc_Class_",cl,".tif",sep=""))
  
}


# ---------------------------------------------------------------------------------- #


# Read all files
# DON'T FORGET TO CHANGE THE DIRECTORY
rstList <- list.files("C:/Users/JG/Desktop/maps", pattern=".tif$", full.names = TRUE)
rstList <- rstList[-1] # Exclude the no-transition class
rstList <- rstList[c(5:12,1:4)] # Re-order the vector with raster file names



# Represent change/no-change maps
rstChangeOnly <- rast(rstList[1]) 

plot(rstChangeOnly) # No change plot

plot(100 - rstChangeOnly) # Reversed/ changed areas %



rstStack <- rast(rstList) # Create a raster stack with all raster percents

# Names for each transition
names(rstStack) <- c(
                      "Deforestation",
                      "Afforestation",
                      "Reforestation",
                      "Forest plantation expansion",
                      "Native vegetation interchange",
                      "Biomass gain",
                      "Rocky exposure / Biomass loss",
                      "Renaturalization",
                      "Farming expansion",
                      "Urban-mining expansion",
                      "Water increase",
                      "Water decrease")


plot(rstStack)


# Validation sum (ignore!!)
rr <- stack(rstList)
rrsum <- calc(rr, sum)
zz <- na.omit(values(rrsum))
unique(zz)
plot(rrsum)




rasterVis::levelplot(rstStack, par.settings=BuRdTheme())

rasterVis::levelplot(rstStack[[9]], par.settings=())


