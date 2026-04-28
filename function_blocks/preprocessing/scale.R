program_block_PP <- function(data, path_og_dataset='', omic='') {

  scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }

  if (omic == 'ref_scRNA') {
    data = lapply(data, function(x) {
      list(counts=scale_matrix(x$counts),
           metadata=x$metadata)})
  } else {
    data = scale_matrix(data)
  }
   
  return(data) 
}
