program_block_DE <- function(uni_data,path_og_dataset='') {
  
  remotes::install_github('xuranw/MuSiC')
  library(MuSiC)
  library(SingleCellExperiment)
  
  bulk.mtx = uni_data$mix 
  
    #convert into single cell experiment format
    sce_object <- SingleCellExperiment(
      assays = list(counts = uni_data$ref_scRNA$ref_sc_peng$counts),
      colData = DataFrame(
        cell_type = uni_data$ref_scRNA$ref_sc_peng$metadata$cell_type,  # Matches clusters argument
        sample = uni_data$ref_scRNA$ref_sc_peng$metadata$sample     # Matches samples argument
      )
    )
    
    uni_pred = music_prop(bulk.mtx = bulk.mtx, sc.sce = sce_object, clusters = 'cell_type',
                          samples = 'sample', select.ct = unique(sce_object$cell_type), verbose = F)
    
    return(t(uni_pred$Est.prop.weighted))
    
}

