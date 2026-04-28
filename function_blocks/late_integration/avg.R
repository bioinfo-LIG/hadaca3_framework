program_block_li <- function(prop1,prop2,path_dataset) { 

  prop = Reduce(`+`, list(prop1,prop2)) / length(list(prop1,prop2))

  return(prop)
}
