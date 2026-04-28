##################################################################################################
### PLEASE only edit the program function between YOUR CODE BEGINS/ENDS HERE                   ###
##################################################################################################

#' The function to predict the A matrix
#' In the provided example, we use basic non-negative least squares (package "nnls"), which consists in minimizing the error term $||Mix - Ref \times Prop||^2$ with only positive entries in the prop matrix.
#'
#' @param mix a matrix of bulk samples (columns) and features (rows)
#' @param ref a matrix of pure cell types (columns) and features (rows)
#' @param ... other parameters that will be ignored
#' 
#' @return the predicted A matrix
#' @examples
#' 
program = function(mix=NULL, ref=NULL, ...) {
  ##
  ## YOUR CODE BEGINS HERE
  ##

  #### RNA-SEQ DATA
  
  # Sequence depth normalization
  seq_depth_normalization <- function(mat) {
    sweep(mat, 2, colSums(mat), "/") * 10^6
  }
  
  mix_normed <- seq_depth_normalization(mix)
  ref_normed <- seq_depth_normalization(ref)
  
  
  # Select only top variable genes
  determine_variable_genes <- function(cpm, n_genes)  {
    
    constant <- median(cpm[cpm > 0])
    vst <- log2(cpm + constant)
    
    gene_variance <- apply(vst, 1, var)
    var_sorted_genes <- names(sort(gene_variance, decreasing = TRUE))
    
    subset_genes <- var_sorted_genes[1:n_genes]
    
    return(subset_genes)
  }
  
  top_genes <- determine_variable_genes(ref_normed, 5000)
  subset_mix <- mix_normed[top_genes, ]
  subset_ref <- ref_normed[top_genes, ]
  
  
  # Now scale the rows so ref and mix are on the same scale
  row_mean <- rowMeans(subset_mix) + rowMeans(subset_ref) + 1e-4
  subset_mix_rna <- sweep(subset_mix, 1, row_mean, '/')
  subset_ref_rna <- sweep(subset_ref, 1, row_mean, '/')
  
  mix <- subset_mix_rna
  ref <- subset_ref_rna
  
  # Perform deconvolution
  rna_prop <- apply(mix, 2, function(b, A) {
    tmp_prop = nnls::nnls(b = b, A = A)$x
    tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
    return(tmp_prop)
  }, A = ref)
  
  rownames(rna_prop) <- colnames(ref)
  
  
  #### METHYLATION DATA
  
  # Filter sites by variance
  select_top_variance_sites <- function(met, threshold) {
    
    variance <- apply(met, 1, var)
    
    var_sorted_sites <- names(sort(variance, decreasing = TRUE))
    # high_var_sites <- var_sorted_sites[1:n_sites]
    high_var_sites <- names(which(variance > threshold))
    
    return(high_var_sites)
  }
  
  high_var_sites <- select_top_variance_sites(ref_met, 0.1)
  
  filtered_mix_met <- mix_met[high_var_sites, ]
  filtered_ref_met <- ref_met[high_var_sites, ]
  
  mix <- filtered_mix_met
  ref <- filtered_ref_met
  
  # Perform deconvolution
  met_prop <- apply(mix, 2, function(b, A) {
    tmp_prop = nnls::nnls(b = b, A = A)$x
    tmp_prop = tmp_prop / sum(tmp_prop) # Sum To One
    return(tmp_prop)
  }, A = ref)
  
  rownames(met_prop) <- colnames(ref)
  
  # Late integration
  prop <- (met_prop + rna_prop) / 2
  
  bad_rna_cell_types <- c("classic", "basal")
  prop[bad_rna_cell_types, ] <- met_prop[bad_rna_cell_types, ]
  prop <- sweep(prop, 2, colSums(prop), "/")  # need to apply scaling per column again

  return(prop)
  
  ##
  ## YOUR CODE ENDS HERE
  ##
}