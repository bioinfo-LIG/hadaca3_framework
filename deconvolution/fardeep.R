

program_block_DE <- function(uni_data,path_og_dataset='') {
  
  if ( !( "FARDEEP" %in% installed.packages() ) ) {
    install.packages("FARDEEP", repos = "http://cran.us.r-project.org", dependencies = TRUE, INSTALL_opts = '--no-lock')
  }
  library(FARDEEP)
  
  
  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref = uni_data$ref[idx_feat,]
  uni_pred = t(FARDEEP::fardeep(
    X = uni_data$ref,
    Y = uni_data$mix)$relative.beta
  )
  
  return(uni_pred) 
}

