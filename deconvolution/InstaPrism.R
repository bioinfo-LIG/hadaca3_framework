

program_block_DE <- function(uni_data,path_og_dataset='') {
  
  library(BayesPrism)
  library(pbapply)
  library(dplyr)
  library(NMF)
  
  idx_feat = intersect(rownames(uni_data$mix), rownames(uni_data$ref))
  uni_data$mix = uni_data$mix[idx_feat,]
  uni_data$ref = uni_data$ref[idx_feat,]
  
  
  state_labels <- colnames(uni_data$ref)
  type_labels <- state_labels
  
  type_labels[grepl("classic", type_labels)] <- "tumor"
  type_labels[grepl("basal", type_labels)] <- "tumor"
  
  
  prism <- BayesPrism::new.prism(
    reference = t(uni_data$ref),
    mixture = base::t(uni_data$mix),
    input.type = "GEP",
    cell.state.labels = state_labels,
    cell.type.labels = type_labels,
    key = "tumor") # create prism obj
  
  
  
  InstaPrism_obj = InstaPrism::InstaPrism(prismObj = prism, input_type = "prism")
  uni_pred = InstaPrism_obj@Post.ini.cs@theta
  
  return(uni_pred) 
}
















prism.states <- function(dat, ref_profiles, nCores = 32) {
  ncores <- nCores - 1
  # define types and states
  state_labels <- colnames(ref_profiles)
  type_labels <- state_labels
  
  ## Define variable types (tumor types in our case)
  ## BrCL1: types and states are equals, tumoral type label is "tumor"
  ## PaCL1: 2 tumoral states, 1 tumoral type
  type_labels[grepl("TUM_", type_labels)] <- "tumor"
  ## PaCL2 and PaPB: 2 tumoral states, 1 tumoral type
  type_labels[grepl("Cancer", type_labels)] <- "tumor"
  ## LuCL: 1 tumoral state and type
  type_labels[grepl("A549", type_labels)] <- "tumor"
  ## BrCL2: 3 tumoral states, 1 tumoral type
  type_labels[grepl("BT474", type_labels)] <- "tumor"
  type_labels[grepl("MCF7", type_labels)] <- "tumor"
  type_labels[grepl("T47D", type_labels)] <- "tumor"
  
  prism <- BayesPrism::new.prism(
    reference = base::t(ref_profiles),
    mixture = base::t(dat),
    input.type = "GEP",
    cell.state.labels = state_labels,
    cell.type.labels = type_labels,
    key = "tumor") # create prism obj
  
  res <- InstaPrism(prismObj = prism, input_type = "prism",
                    n.core = ncores) # run deconv
  A_state <- t(res@Post.ini.cs@theta) # get state props before update
  A_type <- t(res@Post.ini.ct@theta) # get type props after update
  # get results
  ## we need to "deaggregate" tumoral type to tumoral states
  ### first get tumoral states proportions within the tumoral type
  
  ## Deaggregate variable types
  state_labels <- base::colnames(A_state)
  tumoral_states_mask <- state_labels == "tumor" |
    grepl("TUM_", state_labels) |
    grepl("Cancer", state_labels) |
    grepl("A549", state_labels) |
    grepl("BT474", state_labels) |
    grepl("MCF7", state_labels) |
    grepl("T47D", state_labels)
  tumoral_states_labels <- state_labels[tumoral_states_mask] # use of labels rather than mask to make sure order in the dimensions name is not important
  not_tumoral_states_labels <- state_labels[!tumoral_states_mask]
  tumoral_states_prop <- A_state[, tumoral_states_labels]
  if (!is.array(tumoral_states_prop)) {
    tumoral_states_prop <- matrix(tumoral_states_prop, ncol = 1,
                                  dimnames = list(base::rownames(A_state), tumoral_states_labels))
  }
  tumoral_states_prop <- tumoral_states_prop / apply(tumoral_states_prop, 1, sum)
  #### then dispatch tumoral type final proportion to states accordingly
  type_labels <- base::colnames(A_type)
  tumoral_type_mask <- type_labels == "tumor"
  tumoral_type_labels <- type_labels[tumoral_type_mask]
  not_tumoral_type_labels <- type_labels[!tumoral_type_mask]
  A_matrix <- array(dim = dim(A_state), dimnames = dimnames(A_state))
  A_matrix[, tumoral_states_labels] <- A_type[, tumoral_type_labels] * tumoral_states_prop
  ### other states remain unchanged
  A_matrix[, not_tumoral_states_labels] <- A_type[, not_tumoral_type_labels]
  return(t(A_matrix[, colnames(ref_profiles)]))
}