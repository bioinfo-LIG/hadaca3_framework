## Authors: Francisco Avila Cobos (https://github.com/favilaco/deconv_benchmark)
## Adapted by magali richard
## magali.richard@univ-grenoble-alpes.fr
##
## ------------------------------------

generator_hadaca3 <- function(sc_data, phenoData, proportion, Num.mixtures = 1000, pool.size = 100, min.percentage = 1, max.percentage = 99, seed = 24, type = c("rna","met")){ 
  
  CT = unique(phenoData$cellType)
  ?stopifnot(length(CT) >= 2)
  
  set.seed(seed)
  require(dplyr)
  require(gtools)
  
  cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE) #nb of cells type per CT
  colnames(cell.distribution) = c("CT","max.n") # number of cell per cell type is the maximum number of cell we will be able to pick from later on 
  
  Tissues = list()
  Proportions = list()
  
  for(y in 1:Num.mixtures){ 
    
    p = proportion[,y] # get proportion for the given sample
    P = data.frame(CT = names(p), expected = as.numeric(p), stringsAsFactors = FALSE) # get to the right formation for the rest of the function 
    
    
    # Using info in P to build T simultaneously
    chosen_cells <- sapply(which(P$expected != 0), function(x){ 
      
      n.cells = P$expected[x] * pool.size # from proportion get number of cell? 
      #chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]], n.cells)
      chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]], n.cells, replace = TRUE) #mag change
      
      chosen
    }) %>% unlist() # generate a list of "pool size" cells with proportion P according to cell type. 
    
    if (type =="rna"){
      T <- Matrix::rowSums(sc_data[,colnames(sc_data) %in% chosen_cells]) %>% as.data.frame() # somme des expressions des cellules selectionnées pour chaque gene dans le ds single cell.
    }
    if (type =="met"){
      T <- Matrix::rowMeans(sc_data[,colnames(sc_data) %in% chosen_cells]) %>% as.data.frame() # somme des expressions des cellules selectionnées pour chaque gene dans le ds single cell. 
    }
    colnames(T) = paste("mix",y,sep="")
    
    Tissues[[y]] <- T
    
    rm(list=c("T","P","chosen_cells"))
    
  }
  
  T = do.call(cbind.data.frame, Tissues)
  
  
  return(T)
} 

