

fl <- list.files("C:/Users/JG/Google Drive (up200103770@g.uporto.pt)/CIBIO_Desktop/R_projects/TeseMariana_CIBIOPC/OUT/MODS_CIBIOPC/R01/_TRAIN_DATASETS",
           pattern="_SelectedVars_v1.rds$", full.names=TRUE, recursive = TRUE)

selVars <- list()

for(fn in fl){
  
  print(fn)
  fname <- gsub("_SelectedVars_v1.rds","",basename(fn))
  
   tmp <- readRDS(fn)
   
   attr(tmp,"vimpDF") <- NULL
   selVars[[fname]] <-tmp
}

lapply(selVars, FUN = function(x) paste(x,collapse=", "))

speciesVars <- 
data.frame(
speciesNames = names(selVars),
varList = unlist(lapply(selVars, FUN = function(x) paste(x,collapse=", "))))


write.csv(speciesVars, "C:/Users/JG/Desktop/varList_v1.csv", row.names = FALSE)




