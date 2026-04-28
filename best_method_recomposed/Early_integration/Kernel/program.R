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
  #lognorm
  mix_rna = exp(Seurat::LogNormalize(mix_rna))
  ref_bulkRNA = exp(Seurat::LogNormalize(ref_bulkRNA))


  #### FS RNA
  #Toastbulknbfs
  nb_fs_rna = 1e3
  idx_feat = intersect(rownames(mix_rna), rownames(ref_bulkRNA))
  mix_rna = mix_rna[idx_feat,]
  ref_bulkRNA = ref_bulkRNA[idx_feat,]
  hvg <- TOAST::findRefinx(ref_bulkRNA, nmarker = min(nrow(ref_bulkRNA), nb_fs_rna))
  mix_rna = mix_rna[hvg,]
  ref_bulkRNA = ref_bulkRNA[hvg,]


  ###################
  ##### preprocess MET
  #Scale
  scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }
  mix_met = scale_matrix(mix_met)
  ref_met = scale_matrix(ref_met)

  ##### fs MET   
  # SPLSDA
  
  highly_var_idx <- TOAST::findRefinx(ref_met, nmarker=1e4)
  ref_met = ref_met[highly_var_idx,]
  sda <- mixOmics::splsda(t(log(ref_met/(1-ref_met))), Y=colnames(ref_met),
                    keepX = c(1000,1000))
  choose_markers_met <- unique(c(mixOmics::selectVar(sda, comp = 1)$name,
                                  mixOmics::selectVar(sda, comp = 2)$name))

  ref_met = ref_met[choose_markers_met,]
  mix_met = mix_met[choose_markers_met,]



  ######################
  ### Early integration
  #  Kernel


  if (!("mixKernel" %in% installed.packages())) {
    BiocManager::install("mixKernel")
  }
  library(mixKernel)

  # add columns' names and order in the same way
  if (is.null(colnames(mix_rna))) {
    colnames(mix_rna) = paste0("Sample",seq(ncol(mix_rna)))
  }
  if (is.null(colnames(mix_met))) {
    colnames(mix_met) = paste0("Sample",seq(ncol(mix_met)))
  }
  mix_rna = mix_rna[,colnames(mix_met)]


  common_genes <- intersect(rownames(mix_rna), rownames(ref_bulkRNA))

  mix_rna_aligned <- mix_rna[common_genes, ]
  ref_rna_aligned <- ref_bulkRNA[common_genes, ]
  
  common_meth_probes = intersect(rownames(mix_met), rownames(ref_met))
  
  mix_met_aligned <- mix_met[common_meth_probes, ]
  ref_met_aligned <- ref_met[common_meth_probes, ]


  # compute kernels
  kernel_rna <- compute.kernel(t(as.matrix(cbind(mix_rna_aligned,ref_rna_aligned))), kernel.func = "abundance")
  kernel_met <- compute.kernel(t(as.matrix(cbind(mix_met_aligned,ref_met_aligned))), kernel.func = "abundance")
  kernel_all <- combine.kernels(kernel_rna = kernel_rna, kernel_met = kernel_met)
  kernel_pca <- kernel.pca(kernel_all, ncomp = 10)
  
  # Latent space is in kernel_pca$variates$X
  projection = t(kernel_pca$variates$X)
  projection = projection - min(projection)
  mix = projection[,colnames(mix_rna)]
  ref = projection[,colnames(ref_bulkRNA)]

  ##### decovolution    
  ##epic 
  idx_feat = intersect(rownames(mix), rownames(ref))
  
  prop <- t(EPIC::EPIC(as.matrix(mix),
        reference = list(refProfiles = ref, sigGenes = rownames(mix)),
        scaleExprs = F, withOtherCells = F)$mRNAProportions)

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
