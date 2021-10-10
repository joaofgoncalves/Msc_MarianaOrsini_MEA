

ageoffile <- function(x){
  dt = difftime(Sys.time(), file.info(x)$mtime, units="mins")
  return(dt)
}

for(i in 1:10000){
  
  fl <- list.files("C:/temp/RtmpiA7SMY", full.names = TRUE, pattern = ".tif$")
  
  for(fn in fl){
    
    fage <- ageoffile(fn)
    print(fage)
    
    if(fage >= 4){
      cat("Removing file:\n")
      cat(paste(" ->",fn,"\n\n"))
      file.remove(fn)
    }
  }
  cat("Sleeping now....\n\n\n")
  Sys.sleep(60)
}
