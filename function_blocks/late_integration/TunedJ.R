program_block_li <- function(prop1,prop2,path_dataset) { 
  
  prop <- Reduce(`+`, list(prop1,prop2)) / length(list(prop1,prop2))
  bad_rna_cell_types <- c("classic", "basal")
  prop[bad_rna_cell_types, ] <- prop2[bad_rna_cell_types, ]
  prop <- sweep(prop, 2, colSums(prop), "/")  # STO

  
  return(prop)
}
