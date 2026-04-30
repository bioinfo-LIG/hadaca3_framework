# ppID & scpseudobulk & sccluster & SCcluster & ppID & mostmethylated & concatnoscale & RLRpoisson


program <- function(mix_rna=NULL, ref_bulkRNA=NULL, 
                    mix_met=NULL, ref_met=NULL, ref_scRNA=NULL) {


  # Library
  install.packages(c("dplyr","MASS","caret"))
  

  ############### SC RNA 
  #### pp  Preprocess + features selection scRNA. 
  metadata <- do.call(rbind, lapply(ref_scRNA, function(x) x$metadata))
  ref_scRNA <- lapply(ref_scRNA, function(x) as.matrix(x$counts))
  shared_genes <- Reduce(intersect, lapply(ref_scRNA, rownames))
  ref_scRNA <- lapply(ref_scRNA, function(x) x[shared_genes,])
  ref_scRNA <- lapply(ref_scRNA, function(x) as(x,'dgCMatrix'))
  ref_scRNA_all <- do.call(cbind, ref_scRNA)
  counts = ref_scRNA_all
  
  metadata$dataset = sapply(rownames(metadata), function(x) strsplit(x, ".", fixed=T)[[1]][1])
  rownames(metadata) <- sub("^[^.]+\\.", "", rownames(metadata))

  ######################## RNA
  #### PP RNA 
  # Log Norm
  mix_rna = exp(Seurat::LogNormalize(mix_rna))
  ref_bulkRNA = exp(Seurat::LogNormalize(ref_bulkRNA))

  #### FS RNA
  #Scpseudobulk
  t_statistics <- lapply(combn(unique(metadata$cell_type), 2, simplify = F), function(cts) {
      lapply(unique(metadata$dataset), function(ds) {
        idx_col = metadata$dataset == ds & metadata$cell_type %in% cts
        # if (!any(idx_col) ) {
        #   # message(cts," ", ds ,"NULL")
        #   return(NULL)
        # }
        counts_trunc = counts[, which(idx_col)]

      if (is.null(dim(counts_trunc)) || any(dim(counts_trunc) == 0)) {
        return(NULL)
      }

        nb_celltype = length(unique(metadata[colnames(counts_trunc), "cell_type"]))

        apply(counts_trunc, 1, function(x) {

          if (nb_celltype > 1 && length(unique(x)) > 0) {
              res <- try(t.test(x ~ metadata[colnames(counts_trunc), "cell_type"])$statistic, silent = T)
            if (class(res) == "try-error") {return(0)}
            return(res)
          } else {return(NULL)}
        })
      })
    })


  n_top_genes = 20
  top_genes <- unique(do.call(c, lapply(t_statistics, function(x) {
    x <- x[!sapply(x, is.null)]
    Reduce(intersect, lapply(x, function(y) {
      y <- y[is.finite(y)]
      y <- y[order(y)]
      c(names(head(y, n_top_genes)), names(tail(y, n_top_genes)))}))
  })))
  
  res <- aggregate(as.data.frame(t(as.matrix(counts[top_genes, ]))),
                   list(metadata[, "cell_type"]),
                   mean)
  rownames(res) <- res[, 1]
  res <- t(res[-1])
  
  # ref_scRNA = lapply(ref_scRNA_all, function(x) list(counts = x$counts[rownames(res),], metadata = x$metadata))
  # which(rownames(res) %in% rownames(data))
  mix_rna = mix_rna[which(rownames(res) %in% rownames(mix_rna)),]
  ref_bulkRNA = ref_bulkRNA[which(rownames(res) %in% rownames(ref_bulkRNA)),]

   

    # data = data[rownames(res),]



  

  ###################
  ##### preprocess MET
  #ID

  ##### fs MET   
  # Toastpercent

  TOAST_percent_met = 0.8


  nmarker_percent = round(TOAST_percent_met*nrow(ref_met))
  hvp <- TOAST::findRefinx(ref_met, nmarker = nmarker_percent)


  mix_met = mix_met[hvp,]
  ref_met = ref_met[hvp,]




  ######################
  ### Early integration
  # concatscale
  mix = rbind(mix_rna, mix_met)
  ref = rbind(ref_bulkRNA, ref_met)
  
  mix = apply(mix, 2, function(x) {
    tr1 = x - mean(x, na.rm=T)
    pnorm(tr1/sd(x, na.rm=T))
  })
  ref = apply(ref, 2, function(x) {
    tr1 = x - mean(x, na.rm=T)
    pnorm(tr1/sd(x, na.rm=T))
  })

  ##### decovolution    
  ## RLR 
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
