program_block_PP <- function(data, path_og_dataset='', omic='') {
  
  # Sequence depth normalization for RNA
  seq_depth_normalization <- function(mat) {
    sweep(mat, 2, colSums(mat), "/") * 10^6
  }
  
  if (omic == 'ref_scRNA' ) {
    data <- lapply(data, function(x)
      list(counts=seq_depth_normalization(x$counts), metadata=x$metadata))
  } else {
    data <- seq_depth_normalization(data)
  }

  
  return(data) 
}