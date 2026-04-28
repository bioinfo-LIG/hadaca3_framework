# ppID & scpseudobulk & sccluster & SCcluster & ppID & mostmethylated & concatnoscale & RLRpoisson


program <- function(mix_rna=NULL, ref_bulkRNA=NULL, 
                    mix_met=NULL, ref_met=NULL, ref_scRNA=NULL) {


  # Library
  # install.packages(c("dplyr","MASS","caret"))
  # install.packages("TOAST", repos = "http://cran.us.r-project.org", dependencies = TRUE, INSTALL_opts = "--no-lock")
  

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
  mix_met = clean_matrix(mix_met)
  ref_bulkRNA <- clean_matrix(ref_bulkRNA)
  ref_met <- clean_matrix(ref_met)
  ref_scRNA <- lapply(ref_scRNA, function(x) {
        list(counts=as(clean_matrix(x$counts),'dgCMatrix'),
              metadata=x$metadata)}
        )   


  ####################
  ##### prepocess scRNA :
  # LogNorm 
  # ref_scRNA <- lapply(ref_scRNA, function(x) list(counts=Seurat::LogNormalize(x$counts), metadata=x$metadata))

  ##### feature selection scRNA : 
  # ID


  #max_discriminent
  library(dplyr) 
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

  # ref_bulkRNA = ref_bulkRNA[intersect(rownames(ref_bulkRNA),
  #                                           rownames(data)),]
  # data = data[rownames(ref_bulkRNA),]
  hvg <- TOAST::findRefinx(ref_bulkRNA, nmarker = min(nrow(ref_bulkRNA), nb_fs_rna))
  mix_rna = mix_rna[hvg,]
  ref_bulkRNA = ref_bulkRNA[hvg,]
  

  ##### decovolution RNA 
  idx_feat = intersect(rownames(mix_rna), rownames(ref_bulkRNA))
  mix_rna = mix_rna[idx_feat,]
  ref_bulkRNA = ref_bulkRNA[idx_feat,]

  get_weights_poisson <- function(ref) {
    weights = 1/apply(ref, 1, mean)
    return(weights/sum(weights))
  }

  compute_rlr_weighted <- function(beta.m, ref.m, weights, correlation_cutoff = 0.99) {
  library(MASS)
  library(caret)
    est.m <- matrix(nrow = ncol(beta.m), ncol = ncol(ref.m))
    colnames(est.m) <- colnames(ref.m)
    rownames(est.m) <- colnames(beta.m)

    for (s in seq_len(ncol(beta.m))) {
      current_ref <- ref.m
      ref_cols <- colnames(ref.m)
      success <- FALSE

      while (!success) {
        # Try fitting the robust linear model
        tryCatch({
          rlm.o <- MASS::rlm(beta.m[, s] ~ current_ref - 1, maxit = 50, weights = weights)
          coef.v <- summary(rlm.o)$coef[, 1]

          # Normalize: remove negatives and enforce sum to 1
          coef.v[coef.v < 0] <- 0
          total <- sum(coef.v)
          coef.v <- coef.v / total

          # Store results
          est.m[s, ] <- NA
          est.m[s, colnames(current_ref)] <- coef.v
          success <- TRUE
        }, error = function(e) {
          message(paste("RLM failed for sample", s, "- removing correlated columns..."))

          # Remove highly correlated columns
          cor_mat <- cor(current_ref)
          to_remove <- caret::findCorrelation(cor_mat, cutoff = correlation_cutoff)

          if (length(to_remove) == 0) {
            stop("RLM failed and no correlated columns left to remove.")
          }

          current_ref <- current_ref[, -to_remove, drop = FALSE]
        })
      }
    }
    return(t(est.m))
  } 

  prop_RNA = compute_rlr_weighted(beta.m = mix_rna, ref.m = ref_bulkRNA, weights = get_weights_poisson(ref_bulkRNA))
  rownames(prop_RNA) <- colnames(ref_bulkRNA)


  

  ###################
  ##### preprocess MET
  #ID

  ##### fs MET 
  # Maxdiscriminant
  og_ref_met = ref_met
  
  ## 1. For each probe, compute the absolute difference between the probe's methylation level in a cell type versus the global average in the matrix
  ## 2. Ranking probes for each cell type by the score computed above
  ## 3. Extract the top n features per cell type
  n_discriminant <- function(mat, n) {
    mean_ct <- rowMeans(mat)
    df_mat <- mat %>% apply(2, function(x) round((abs(x - mean_ct)), 4))
    colnames(df_mat) %>%
      lapply(function(x) {
        df_mat[, x] %>%
          sort(decreasing = TRUE) %>%
          head(n) %>% names()}) %>% Reduce(c, .)
  }
  ## Return the maximal nb of probes (up to n_max) s.t. there is no overlap across cell types
  least_n_discriminant <- function(mat, n_max) {
    nb_features <- c(1:n_max) %>% parallel::mclapply(function(i) {
      top_n_each <- n_discriminant(mat, i)
      return(max(table(top_n_each)))
    }, mc.cores = 12)
    least_n <- min(which(unlist(nb_features) > 1)) - 1
    return(unique(n_discriminant(mat, least_n)))
  }
  
  features_met <- least_n_discriminant(og_ref_met, 100)
  
  mix_met = mix_met[features_met,]
  ref_met = ref_met[features_met,]


  ##### decovolution MET 
  # function defined in RNA decovolution
  prop_MET = compute_rlr_weighted(beta.m = mix_met, ref.m = ref_met, weights = get_weights_poisson(ref_met))
  rownames(prop_MET) <- colnames(ref_met)




  ###################
  ##### Late integration  
  idx_feat_rna = intersect(rownames(mix_rna), rownames(ref_bulkRNA))
  mix_rna = mix_rna[idx_feat_rna,]
  ref_bulkRNA = ref_bulkRNA[idx_feat_rna,]
  idx_feat_met = intersect(rownames(mix_met), rownames(ref_met))
  mix_met = mix_met[idx_feat_met,]
  ref_met = ref_met[idx_feat_met,]
  

  real_rna = mix_rna / colSums(mix_rna)
  reconstructed_rna = (ref_bulkRNA %*% prop_RNA) / colSums(ref_bulkRNA %*% prop_RNA)
  rmse_rna = sqrt(mean(((as.matrix(reconstructed_rna)) - as.matrix(real_rna))^2))
    
  real_met = mix_met / colSums(mix_met)
  reconstructed_met = (ref_met %*% prop_MET) / colSums(ref_met %*% prop_MET)
  rmse_met = sqrt(mean(((as.matrix(reconstructed_met)) - as.matrix(real_met))^2))
  
  # Normalize RMSEs
  rmse_rna_norm = rmse_rna / (rmse_rna + rmse_met)
  rmse_met_norm = rmse_met / (rmse_rna + rmse_met)
  
  # RMSE as weights  
  prop <- as.matrix( rmse_met_norm * prop_RNA + rmse_rna_norm  * prop_MET )

  
  return(prop)


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
