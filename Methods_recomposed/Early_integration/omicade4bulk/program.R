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
  # SCALE

  scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }

  mix_rna = scale_matrix(mix_rna)
  ref_bulkRNA = scale_matrix(ref_bulkRNA)

  #### FS RNA
  #Toastvst
  
  nb_fs_rna = 1e3

  determine_variable_genes <- function(mat, n_genes=nb_fs_rna) {
    vst <- log2(mat + median(mat[mat > 0]))
    var_sorted_genes <- TOAST::findRefinx(vst, nmarker = min(nrow(ref_bulkRNA), n_genes))
    return(var_sorted_genes)
  }
  top_genes <- determine_variable_genes(ref_bulkRNA, nb_fs_rna)
  top_gene_names <- rownames(ref_bulkRNA)[top_genes]


  mix_rna = mix_rna[ intersect(top_gene_names,rownames(mix_rna)),]
  ref_bulkRNA = ref_bulkRNA[ intersect(top_gene_names,rownames(ref_bulkRNA)),]


  ###################
  ##### preprocess MET
  #Scale
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
  #  omicade4bulk

  if (!("omicade4" %in% installed.packages())) {
    BiocManager::install("omicade4")
  }
  library(omicade4)
  
  # add samples' names
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

  
  # run omicade4
  data_combi <- list(RNA=as.matrix(cbind(mix_rna_aligned,ref_rna_aligned)),
                      DNAm=as.matrix(cbind(mix_met_aligned,ref_met_aligned)))
  combi <- mcia(data_combi, cia.nf=10)
  
  # Latent space is in combi$mcoa$SynVar
  projection = t(combi$mcoa$SynVar)
  projection = projection - min(projection)
  colnames(projection) = colnames(cbind(mix_rna_aligned,ref_rna_aligned))
  
  # Retrieve D and T
  mix = projection[,colnames(mix_rna_aligned)]
  ref = projection[,colnames(ref_rna_aligned)]


  ##### decovolution    
  ##LM 
  idx_feat = intersect(rownames(mix), rownames(ref))
  
  # Estimation of proportions
  prop = apply(mix[idx_feat,], 2, function(b, A) {
    tmp_prop = lm(b ~ A - 1)$coefficients  # Using `-1` to remove the intercept
    tmp_prop[tmp_prop < 0] = 0
    tmp_prop = tmp_prop / sum(tmp_prop)    # Sum To One
    return(tmp_prop)
  }, A=ref[idx_feat,])

  # Labeling of estimated proportions 
  rownames(prop) = colnames(ref)
  

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
