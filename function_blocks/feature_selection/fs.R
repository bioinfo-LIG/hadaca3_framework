program_block_FS <- function(data, path_og_dataset='') {
  
  drop_null_ref_cols <- function(ref_matrix){
    non_null_rows = apply(ref_matrix != 0,MARGIN = 1, FUN = any, simplify = TRUE)
    return(ref_matrix[non_null_rows,])
   }

  
  data <- drop_null_ref_cols(data)

  return(data) 
}

