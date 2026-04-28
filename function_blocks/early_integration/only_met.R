program_block_EI <- function(rna_unit,met_unit,path_dataset) { 

  within(met_unit, ref_scRNA <- rna_unit$ref_scRNA)
  return(met_unit)
}
