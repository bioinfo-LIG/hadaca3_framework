# ppID & scpseudobulk & sccluster & SCcluster & ppID & mostmethylated & concatnoscale & RLRpoisson


program <- function(mix_rna=NULL, ref_bulkRNA=NULL, 
                    mix_met=NULL, ref_met=NULL, ref_scRNA=NULL) {

    # preprocess 
    scale_matrix <- function(mat) {
    mat = sweep(mat, 2, colSums(mat), "/")
    return(mat)
  }

  if (omic == 'ref_scRNA') {
    data = lapply(data, function(x) {
      list(counts=scale_matrix(x$counts),
           metadata=x$metadata)})
  } else {
    data = scale_matrix(data)
  }
  # Creation of an index, idx_feat, corresponding to the intersection of features present in the references and those present in the mixtures.
  idx_feat = intersect(rownames(mix), rownames(ref))
  
  # Estimation of proportions
  prop = apply(mix[idx_feat,], 2, function(b, A) {
    tmp_prop = lm(b ~ A - 1)$coefficients  # Using `-1` to remove the intercept
    # tmp_prop = nnls::nnls(b=b,A=A)$x  
    tmp_prop[tmp_prop < 0] = 0
    tmp_prop = tmp_prop / sum(tmp_prop)    # Sum To One
    return(tmp_prop)
  }, A=ref[idx_feat,])

  # Labeling of estimated proportions 
  rownames(prop) = colnames(ref)
  return(prop)
  

}