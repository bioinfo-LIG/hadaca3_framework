# ppID & scpseudobulk & sccluster & SCcluster & ppID & mostmethylated & concatnoscale & RLRpoisson


program <- function(mix_rna=NULL, ref_bulkRNA=NULL, 
                    mix_met=NULL, ref_met=NULL, ref_scRNA=NULL) {


  # Library
  install.packages(c("dplyr","MASS","caret"))
  

  #max_discriminent
  library(dplyr) 
  # Decovolution 
  library(MASS)
  library(caret)

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


  library(Seurat)
    





  ##### preprocess scRNA
  seurat_single = lapply(ref_scRNA, function(x) CreateSeuratObject(x$counts, meta.data=x$metadata))
  
  # Merge all
  seurat_merge = Reduce(merge, seurat_single)
  DefaultAssay(seurat_merge) <- "RNA"
  seurat_merge <- JoinLayers(seurat_merge, overwrite = T)
  seurat_merge@meta.data$dataset = unlist(lapply(seq_along(ref_scRNA), function(x)
    rep(names(ref_scRNA)[x], ncol(ref_scRNA[[x]]$counts))))
  
  options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1 GB
  
  process_scrna = function(obj) {
    obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000, verbose=F) 
    obj <- ScaleData(object = obj, features = rownames(obj@assays$RNA)) 
    obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose=F)
    obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca", verbose=F)
    obj <- FindClusters(obj, resolution = 0.5, cluster.name = "clust_res.5")
    obj@meta.data[['clust_res.5']] = Idents(obj)
    Idents(obj) = obj@meta.data$cell_type[match(Idents(obj), obj@meta.data$clust_res.5)]  
    return(obj)
  }
  seurat_merge = process_scrna(seurat_merge)
  # sc = list("ref_cluster"=list(counts=seurat_merge@assays$RNA$counts,
  #                               metadata=data.frame(cell_type=seurat_merge$cell_type,
  #                                                   sample=seurat_merge$sample,
  #                                                   dataset=seurat_merge$dataset),
  #                               seurat_clustered=seurat_merge) )

  ##### feature selection scRNA : 
    

    sc_markers = FindAllMarkers(seurat_merge, assay = NULL, features = NULL,
                              logfc.threshold = 0.1, test.use = "wilcox", slot = "data")
    sc_markers = sc_markers[which(sc_markers$p_val_adj < 0.05 & sc_markers$pct.1>0.6 & sc_markers$pct.2<0.3), ]$gene


    counts = seurat_merge@assays$RNA$counts[sc_markers,]
    metadata = data.frame(cell_type=seurat_merge$cell_type,
                                                    sample=seurat_merge$sample,
                                                    dataset=seurat_merge$dataset)


  ##### PP RNA
  ## ID

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



  ####################
  ##### prepocess scRNA :
  # ID 


  # ID
  

  ###################
  ##### preprocess MET
  #ID

  ##### fs MET   # most methylated
  
  # highest methylation based on biological knowledge (cancer has high proportion of methylated probes)
  quantiles = apply(ref_met, 2, quantile, 0.75) # select 25% of the probes
  methylated_probes = Reduce(union, lapply(1:length(quantiles), function(x)  
    rownames(ref_met)[which(ref_met[,x] > quantiles[x])]))

  mix_met = mix_met[methylated_probes,]
  ref_met = ref_met[methylated_probes,]




  ######################
  ### Early integration
  # concatnoscale
  mix = rbind(mix_rna, mix_met)
  ref = rbind(ref_bulkRNA, ref_met)


  ##### decovolution  

  idx_feat = intersect(rownames(mix), rownames(ref))
  mix = mix[idx_feat,]
  ref = ref[idx_feat,]

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



  prop = compute_rlr_weighted(beta.m = mix, ref.m = ref, weights = get_weights_poisson(ref))
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
