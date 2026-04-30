# ppID & scpseudobulk & sccluster & SCcluster & ppID & mostmethylated & concatnoscale & RLRpoisson


program <- function(mix_rna=NULL, ref_bulkRNA=NULL, 
                    mix_met=NULL, ref_met=NULL, ref_scRNA=NULL) {


  # Library
  install.packages(c("dplyr","MASS","caret"))
  

  clean_matrix <- function(matrix){
      matrix <- as.matrix(matrix)
      # Remove rows  with NA or Inf
      matrix <- matrix[complete.cases(matrix), , drop = FALSE]
      matrix <- matrix[rowSums(is.infinite(matrix)) == 0, , drop = FALSE]
      # Remove rows with all zeros
      matrix <- matrix[rowSums(matrix) > 0, ]  # Remove rows with all zeros
      # Remove rows with zero variance
      matrix <- matrix[apply(matrix, 1, var) > 0, , drop = FALSE]
      return(matrix)
  }


  mix_rna = clean_matrix(mix_rna)
  ref_bulkRNA <- clean_matrix(ref_bulkRNA)



  ####################
  ##### prepocess scRNA :
  # LogNorm 
  # ref_scRNA <- lapply(ref_scRNA, function(x) list(counts=Seurat::LogNormalize(x$counts), metadata=x$metadata))

  ##### feature selection scRNA : 
  # ID

  # Decovolution 
  library(MASS)
  library(caret)

  ##### preprocess RNA
  #Log_Norm
  mix_rna = exp(Seurat::LogNormalize(mix_rna))
  ref_bulkRNA = exp(Seurat::LogNormalize(ref_bulkRNA))
  

  ##### FS RNA 
  # Toastbulknbfs
  nb_fs_rna = 1e3
  idx_feat = intersect(rownames(mix_rna), rownames(ref_bulkRNA))
  mix_rna = mix_rna[idx_feat,]
  ref_bulkRNA = ref_bulkRNA[idx_feat,]
  hvg <- TOAST::findRefinx(ref_bulkRNA, nmarker = min(nrow(ref_bulkRNA), nb_fs_rna))
  mix_rna = mix_rna[hvg,]
  ref_bulkRNA = ref_bulkRNA[hvg,]
    

  ###################
  ##### preprocess MET
  # no LogNorm 
  mix_met = exp(Seurat::LogNormalize(mix_met))
  ref_met = exp(Seurat::LogNormalize(ref_met))

  ##### fs MET 
  # ID
  





  ###################
  ##### early integration
  

  ###################
  ## decovolution
  # RLR
  idx_feat = intersect(rownames(mix), rownames(ref))
  mix = mix[idx_feat,]
  ref = ref[idx_feat,]

  library(caret)

  
  # Define the fallback function
  get_epidish_with_fallback <- function(beta.m, ref.m, cutoff = 0.99) {
    safe_epidish <- function(beta, ref) {
      tryCatch({
        return(t(EpiDISH::epidish(beta, ref, method = "RPC")$estF))
      }, error = function(e) {
        warning("EpiDISH failed, attempting to remove correlated columns...")

        cor_matrix <- cor(ref)
        to_remove <- caret::findCorrelation(cor_matrix, cutoff = cutoff)
        if (length(to_remove) == 0) {
          warning("Still failing after correlation removal. Returning zero matrix.")
          return(matrix(0, nrow = ncol(ref), ncol = ncol(beta),
                        dimnames = list(colnames(ref), colnames(beta))))
        }

        ref_clean <- ref[, -to_remove, drop = FALSE]

        # Try again with reduced matrix
        tryCatch({
          return(t(EpiDISH::epidish(beta, ref_clean, method = "RPC")$estF))
        }, error = function(e2) {
          warning("EpiDISH still failed after cleaning. Returning zero matrix.")
          return(matrix(0, nrow = ncol(ref), ncol = ncol(beta),
                        dimnames = list(colnames(ref), colnames(beta))))
        })
      })
    }

    return(safe_epidish(beta.m, ref.m))
  }


  prop <- get_epidish_with_fallback(mix,ref)


  rownames(prop) <- colnames(ref)
  
  return(prop_RNA)


}



install.packages = function (pkgs, repos="https://cloud.r-project.org", ...) {
  installed_packages <- installed.packages( )
  for (package in pkgs ) {
    if ( !{ package %in% installed_packages } ) {
     print(x = paste("Installation of ", package, sep = "") )
      utils::install.packages(
        pkgs = package,
        repos = repos,
        ...
      )
    } else {
      print(x = paste(package, " is installed.", sep = "") )
    }
  }
}
