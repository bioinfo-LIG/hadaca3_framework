program_block_DE <- function(uni_data,path_og_dataset='') {
  
  
  library(MuSiC)
  library(SingleCellExperiment)
  
  bulk.mtx = uni_data$mix 
  
  metadata = uni_data$ref_scRNA[[1]]$metadata
  counts = uni_data$ref_scRNA[[1]]$counts

  sce_object <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(
      cell_type = metadata$cell_type,  # Matches clusters argument
      sample = metadata$sample     # Matches samples argument
    )
  )

  uni_pred = music_prop(bulk.mtx = bulk.mtx, sc.sce = sce_object, clusters = 'cell_type',
                        samples = 'sample', select.ct = unique(sce_object$cell_type), verbose = F)
  
  return(t(uni_pred$Est.prop.weighted))
    
}

