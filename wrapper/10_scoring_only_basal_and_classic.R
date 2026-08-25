source("data_processing.R")

weighgeomMean <- function(x,w) {
  prod(mapply(function(a,b) a^b, a=x, b=w), na.rm=T)^(1/sum(w))
}

#########################
# homogenization function to find the best match between real and predicted A matrices for unsupervised methods. -USELESS
homogeneized_cor_mat = function(A_real, A_pred) {
    cmat = cor(t(A_real),t(A_pred))
    pvec <- c(clue::solve_LSAP((1+cmat)^2,maximum = TRUE))
    return(A_pred[pvec,])
}

#########################
# Pre-treatment of predicted A
#not useful because reference based algorithm
prepare_A <- function(A_real, A_pred) {
  N <- ncol(A_real)
  K <- nrow(A_real)
    
  stopifnot(K > 1)
  stopifnot(ncol(A_pred) == N)
  stopifnot(!anyNA(A_pred))
    
  ### STEP 1 : matching the number of predicted components to the real number of cell types K
    
  ## if predicting too few cell types
  if (nrow(A_pred) < K) {
    #set random positive values closed to 0 for missing rows
    set.seed(1)
    random_data = abs(jitter(matrix(data = 0, nrow = K - nrow(A_pred), ncol = N), factor = 0.01))
    A_pred <- rbind(A_pred, random_data)
    print("Add rows of 0s to match K")
  }
           
  ## if predicting too many cell types
  
  ### strategy 1: keep the most abundant cell types
  if (nrow(A_pred) > K) {
    A_pred <- A_pred[order(rowSums(A_pred), decreasing = TRUE)[1:K],]
    print("Number of predicted cell types exceeds K, only the K most abundant cell types are kept for scoring.")
  }
  
  ### strategy 2: clustering of similar cell types to match K components -> SLIM
        
  ### STEP 2 : reordering predicted components to find the best match -> SLIM
  A_pred_best_cor = as.matrix(homogeneized_cor_mat(A_real, A_pred))
  
  return(A_pred_best_cor)
}

#########################
# Global Pearson/Spearman correlation coefficients
correlationP_tot = function(A_real, A_pred) {
  if (all(c(var(A_pred))==0)) {return(-1)} #worst case scenario where all cell types are predicted with the same proportions
  return(cor(c(A_real), c(A_pred), method = "pearson"))
}
correlationS_tot = function(A_real, A_pred) {
  if (all(c(var(A_pred))==0)) {return(-1)} #worst case scenario where all cell types are predicted with the same proportions
  return(cor(c(A_real), c(A_pred), method = "spearman"))
}

#########################
# Mean column/sample Pearson/Spearman correlation coefficients
correlationP_col = function(A_real, A_pred) {
  res = c()
  for (i in seq(ncol(A_real))) {
    sd_pred = sd(A_pred[, i], na.rm = TRUE)
    sd_real = sd(A_real[, i], na.rm = TRUE)
    
    if (!is.na(sd_pred) && !is.na(sd_real) && sd_pred > 0 && sd_real > 0) {
    # if (sd(A_pred[i, ]) > 0 & sd(A_real[i, ]) > 0) {

      res[i] = cor(A_real[, i], A_pred[, i], method = "pearson")
    }
  }
  res = res[!is.na(res)]
  print(paste0(length(res), " columns are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}
correlationS_col = function(A_real, A_pred) {
  res = c()
  for (i in seq(ncol(A_real))) {
    sd_pred = sd(A_pred[, i], na.rm = TRUE)
    sd_real = sd(A_real[, i], na.rm = TRUE)
    
    if (!is.na(sd_pred) && !is.na(sd_real) && sd_pred > 0 && sd_real > 0) {
    # if (sd(A_pred[i, ]) > 0 & sd(A_real[i, ]) > 0) {

      res[i] = cor(A_real[, i], A_pred[, i], method = "spearman")
    }
  }
  res = res[!is.na(res)]
  print(paste0(length(res), " columns are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}

#########################
# Mean row/cell type Pearson/Spearman correlation coefficients
correlationP_row = function (A_real, A_pred) {
  res = c()
  for (i in seq(nrow(A_real))) {
    sd_pred = sd(A_pred[, i], na.rm = TRUE)
    sd_real = sd(A_real[, i], na.rm = TRUE)
    
    if (!is.na(sd_pred) && !is.na(sd_real) && sd_pred > 0 && sd_real > 0) {
    # if (sd(A_pred[i, ]) > 0 & sd(A_real[i, ]) > 0) {
      res[i] = cor(A_real[i, ], A_pred[i, ], method = "pearson")
    }
  }
  res = res[!is.na(res)]
  # print(paste0(length(res), " rows are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}
correlationS_row = function (A_real, A_pred) {
  res = c()
  for (i in seq(nrow(A_real))) {
    sd_pred = sd(A_pred[, i], na.rm = TRUE)
    sd_real = sd(A_real[, i], na.rm = TRUE)
    
    if (!is.na(sd_pred) && !is.na(sd_real) && sd_pred > 0 && sd_real > 0) {
    # if (sd(A_pred[i, ]) > 0 & sd(A_real[i, ]) > 0) {
      res[i] = cor(A_real[i, ], A_pred[i, ], method = "spearman")
    }
  }
  res = res[!is.na(res)]
  # print(paste0(length(res), " rows are kept for correlation analysis"))
  if (length(res)==0) {return(-1)}
  return(mean(res))
}

#########################
# Aitchison distance
eval_Aitchison = function(A_real, A_pred, min = 1e-9) {
  # Aitchison dist doesn't like 0 
  A_real[A_real < min] = min
  A_pred[A_pred < min] = min
  
  # Aitchison apply a transformation to manage compositional data, so only sample by sample comparison have sens
  res = c()
  for (i in seq(ncol(A_real))) {
    res[i] = coda.base::dist(rbind(t(A_real[,i]), t(A_pred[,i])), method = "aitchison")[1]
  }
  res = res[!is.na(res)]
  return(mean(res))
}

#########################
# SDID with a twist : i don't change the sign as described in MPRA paper
# https://mpra.ub.uni-muenchen.de/84387/
eval_SDID = function(A_real, A_pred) {
  cos_angle = sapply(seq(ncol(A_real)), function(x)
    A_real[,x]%*%A_pred[,x]/norm(A_real[,x],type="2")/norm(A_pred[,x],type="2"))
  sin_angle = sapply(cos_angle, function(x) sqrt(1-min(1,x^2)))
  sdid = sqrt(sin_angle)
  return(mean(sdid))
}

#########################
# angular distance
eval_AID = function(A_real, A_pred) {
  angle_radian = sapply(seq(ncol(A_real)), function(x)
    acos(min(1,(A_real[,x]%*%A_pred[,x])/norm(A_real[,x],type="2")/norm(A_pred[,x],type="2"))))
  angle_degree = angle_radian * 90 / (pi/2)
  return(mean(angle_degree))
}

#########################
# Jensen Shannon divergence
eval_JSD = function(A_real, A_pred, min = 1e-9) {
  jsd = sapply(seq(ncol(A_real)), function(x) {
    pred = A_pred[,x]
    pred[pred < min] = min
    n <- 0.5 * (A_real[,x] + pred)
    JS <- 0.5 * (sum(A_real[,x] * log(A_real[,x] / n)) + sum(pred * log(pred / n)))})
  return(mean(jsd))
}

#########################
# RMSE
eval_RMSE = function(A_real, A_pred) {
  return(sqrt(mean((A_real - A_pred)^2)))
}

#########################
# MAE
eval_MAE = function (A_real, A_pred){
  return(mean(abs(A_real - A_pred)))
}

#########################
# Scoring function 
scoring_function <- function(A_real, A_pred) {
  # pre-treatment of predicted A
  #A_pred = prepare_A(A_real = A_real, A_pred = A_pred)
    # If A_pred is all NaN, return all-zero scores
  if (all(is.na(A_pred))) {
    metric_names <- c("score_aggreg",
                      rep(c("pearson_row","spearman_row",
                            "pearson_tot","pearson_col","spearman_tot","spearman_col",
                            "rmse","mae","aitchison","jsd","sdid","aid"), times = 2))
    metric_names[(2+length(c("pearson_row","spearman_row",
                              "pearson_tot","pearson_col","spearman_tot","spearman_col",
                              "rmse","mae","aitchison","jsd","sdid","aid"))):
                   length(metric_names)] <-
      paste0(c("pearson_row","spearman_row",
               "pearson_tot","pearson_col","spearman_tot","spearman_col",
               "rmse","mae","aitchison","jsd","sdid","aid"), '_norm')
    
    return(setNames(rep(0, length(metric_names)), metric_names))
  }
  # scoring with different metrics
  if (nrow(A_pred)==nrow(A_real)) {
    
    # Match cell types
    if (all(rownames(A_pred) %in% rownames(A_real))) {
      A_pred = A_pred[rownames(A_real),]
    }
    else {stop(paste0("Cell types names do not match:\n Predicted cell types are ",paste(rownames(A_pred), collapse=', '),"\n Real cell types are ",paste(rownames(A_real), collapse=', ')))}
    rmse = NA
    mae = NA
    aitchison = NA
    sdid = NA
    aid = NA
    jsd = NA
    pearson_tot = NA
    pearson_col = NA
    pearson_row = correlationP_row(A_real, A_pred[rownames(A_real),])
    spearman_tot = NA
    spearman_col = NA
    spearman_row = correlationS_row(A_real, A_pred[rownames(A_real),])
    weights_spec = c(1/2,1/2,rep(0,10))
  }
  else if (nrow(A_pred)==(nrow(A_real)+1)) { # In case of an extra cell type in the reference matrix compared to the source/ground truth
    # Match cell types
    if (all(rownames(A_real) %in% rownames(A_pred))) {
      extra_type = matrix(0,ncol=ncol(A_real),nrow=1)
      rownames(extra_type) = setdiff(rownames(A_pred),rownames(A_real))
      A_real = rbind(A_real,extra_type)
      A_pred = A_pred[rownames(A_real),]
    }
    else {stop(paste0("Cell types names do not match:\n Predicted cell types are ",paste(rownames(A_pred), collapse=', '),"\n Real cell types are ",paste(rownames(A_real), collapse=', ')))}
    rmse = NA
    mae = NA
    aitchison = NA
    sdid = NA
    aid = NA
    jsd = NA
    pearson_tot = NA
    pearson_col = NA
    pearson_row = correlationP_row(A_real, A_pred[rownames(A_real),])
    spearman_tot = NA
    spearman_col = NA
    spearman_row = correlationS_row(A_real, A_pred[rownames(A_real),])
    weights_spec = c(1/2,1/2,rep(0,10))
  }
  else if (nrow(A_pred)==(nrow(A_real)-1)) { # In case of a missing cell type in the reference matrix compared to the source/ground truth
    # Match cell types
    if (all(rownames(A_pred) %in% rownames(A_real))) {
      A_real = A_real[rownames(A_pred),]
    }
    else {stop(paste0("Cell types names do not match:\n Predicted cell types are ",paste(rownames(A_pred), collapse=', '),"\n Real cell types are ",paste(rownames(A_real), collapse=', ')))}
    rmse = NA
    mae = NA
    aitchison = NA
    sdid = NA
    aid = NA
    jsd = NA
    pearson_tot = NA
    pearson_col = NA
    pearson_row = correlationP_row(A_real, A_pred[rownames(A_real),])
    spearman_tot = NA
    spearman_col = NA
    spearman_row = correlationS_row(A_real, A_pred[rownames(A_real),])
    weights_spec = c(1/2,1/2,rep(0,10))
  }
  else if (nrow(A_pred) > nrow(A_real) & setequal(rownames(A_real), c("basal",'classic'))) { # partial ground truth only for the in vivo dataset
    rmse = NA
    mae = NA
    aitchison = NA
    sdid = NA
    aid = NA
    jsd = NA
    pearson_tot = NA
    pearson_col = NA
    pearson_row = correlationP_row(A_real, A_pred[rownames(A_real),])
    spearman_tot = NA
    spearman_col = NA
    spearman_row = correlationS_row(A_real, A_pred[rownames(A_real),])
    weights_spec = c(1/2,1/2,rep(0,10))
  }
  
else{
    rmse = NA
    mae = NA
    aitchison = NA
    sdid = NA
    aid = NA
    jsd = NA
    pearson_tot = NA
    pearson_col = NA
    pearson_row = 0
    spearman_tot = NA
    spearman_col = NA
    spearman_row = 0
    weights_spec = 0
}
  all_judges = data.frame("pearson_row"=pearson_row,
                          "spearman_row"=spearman_row,
                          "pearson_tot"=pearson_tot,
                          "pearson_col"=pearson_col,
                          "spearman_tot"=spearman_tot,
                          "spearman_col"=spearman_col,
                          "rmse"=rmse,
                          "mae"=mae,
                          "aitchison"=aitchison,
                          "jsd"=jsd,
                          "sdid"=sdid,
                          "aid"=aid)
  
  # generate best/worst possible metrics
  fake_worst_pred = apply(A_real, 2, function(prop) { 
    tmp = rep(1e-9, length(prop))
    tmp[which.min(prop)] = 1
    return(tmp)})
  judge_candidate = rbind(all_judges,
                          data.frame("pearson_row"=1,
                                     "spearman_row"=1,
                                     "pearson_tot"=1,
                                     "pearson_col"=1,
                                     "spearman_tot"=1,
                                     "spearman_col"=1,
                                     "rmse"=0,
                                     "mae"=0,
                                     "aitchison"=0,
                                     "jsd"=0,
                                     "sdid"=0,
                                     "aid"=0),
                          data.frame("pearson_row"=-1,
                                     "spearman_row"=-1,
                                     "pearson_tot"=-1,
                                     "pearson_col"=-1,
                                     "spearman_tot"=-1,
                                     "spearman_col"=-1,
                                     "rmse"=1,
                                     "mae"=1,
                                     "aitchison"=1, #verify it's the correct worst distance
                                     "jsd"=1, #verify it's the correct worst divergence
                                     
                                    #  "rmse"=eval_RMSE(A_real, fake_worst_pred),
                                    #  "mae"=eval_MAE(A_real, fake_worst_pred),
                                    #  "aitchison"=eval_Aitchison(A_real, fake_worst_pred), #verify it's the correct worst distance
                                    #  "jsd"=eval_JSD(A_real, fake_worst_pred), #verify it's the correct worst divergence
                                     "sdid"=1,
                                     "aid"=90 #verify it's the correct worst angle in multi-D
                                     ))
  
  # normalize scores with a linear shift between [0,1]
  # strategy when normalizing only one method
  CenterScaleNorm <- function(x) {
    tr1 = x - min(x, na.rm=T)
    tr2 = tr1/(max(tr1, na.rm=T))
    return(tr2)
  }
  judge_candidate_norm = apply(judge_candidate, 2, CenterScaleNorm)
  
  # transform scores s.t. 1 is the best score
  judge_candidate_norm = 1 - judge_candidate_norm
  judge_candidate_norm[,grep("pearson",colnames(judge_candidate_norm))] = 1 - judge_candidate_norm[,grep("pearson",colnames(judge_candidate_norm))]
  judge_candidate_norm[,grep("spearman",colnames(judge_candidate_norm))] = 1 - judge_candidate_norm[,grep("spearman",colnames(judge_candidate_norm))]

  # Average over judges with the geometric mean for the candidate of interest
  #score_aggreg = exp(mean(log(judge_candidate_norm[1,]),na.rm=T))
  # weights = c(rep(1/4*1/2,2),
  #             rep(1/4*1/4,4),
  #             rep(1/4*1/4,4),
  #             rep(1/4*1/2,2))
  # if (nrow(A_pred) > nrow(A_real) & setequal(rownames(A_real), c("basal",'classic'))) {weights = weights_spec}
  weights = weights_spec
  score_aggreg = weighgeomMean(judge_candidate_norm[1,],
                               weights)
  

# Define the metrics you want
metrics <- c("pearson_row", "spearman_row")

# Create the vector directly
all_scores <- c(
  score_aggreg, 
  pearson_row, spearman_row, 
  judge_candidate_norm[1, metrics]
)

names(all_scores) <- c("score_aggreg", metrics, paste0(metrics, "_norm"))
  return(all_scores)
}




######################################
# read data 
#####################################

files <- system("find .  -maxdepth 1 -name '*.h5'", intern = TRUE)


library(data.table)
library(stringr)

# 1. Keep the original paths untouched
# files <- system(...) 

# 2. Create a "clean" version just for parsing logic
# (This removes './' and '.h5' so we can split by underscores)
fnames_clean <- gsub("^\\./", "", files)
fnames_clean <- gsub(".h5$", "", fnames_clean)

# --- PROCESS LATE INTEGRATION (pred-li) ---
# Find the indices of files that start with pred-li
late_idx <- grep("^pred-li", fnames_clean)

# Get the clean names for splitting
late_names_only <- fnames_clean[late_idx]
late_parts <- str_split_fixed(late_names_only, "_", 15)

late_int <- data.table(
  full_path = files[late_idx], # <--- USES ORIGINAL PATHS
  dataset   = str_extract(late_parts[, 1], "(?<=pred-li-).*"),
  ref       = late_parts[, 2],
  preprocessing_mixRNA = late_parts[, 3],
  feature_selection_mixRNA = late_parts[, 4],
  preprocessing_RNA    = late_parts[, 5],
  feature_selection_RNA    = late_parts[, 6],
  preprocessing_scRNA  = late_parts[, 7],
  feature_selection_scRNA  = late_parts[, 8],
  deconvolution_rna  = late_parts[, 9],
  preprocessing_mixMET = late_parts[, 10],
  feature_selection_mixMET = late_parts[, 11],
  preprocessing_MET    = late_parts[, 12],
  feature_selection_MET    = late_parts[, 13],
  deconvolution_met    = late_parts[, 14],
  late_integration    = late_parts[, 15]
)

# --- PROCESS EARLY INTEGRATION (out-li) ---
# Find the indices of files that start with out-li
early_idx <- grep("^out-li", fnames_clean)

# Get the clean names for splitting
early_names_only <- fnames_clean[early_idx]
early_parts <- str_split_fixed(early_names_only, "_", 14)

early_int <- data.table(
  full_path = files[early_idx], # <--- USES ORIGINAL PATHS
  dataset   = str_extract(early_parts[, 1], "(?<=out-li-).*"),
  ref       = early_parts[, 2],
  preprocessing_mixRNA = early_parts[, 3],
  feature_selection_mixRNA = early_parts[, 4],
  preprocessing_RNA    = early_parts[, 5],
  feature_selection_RNA    = early_parts[, 6],
  preprocessing_scRNA  = early_parts[, 7],
  feature_selection_scRNA  = early_parts[, 8],
  preprocessing_mixMET = early_parts[, 9],
  feature_selection_mixMET = early_parts[, 10],
  preprocessing_MET    = early_parts[, 11],
  feature_selection_MET    = early_parts[, 12],
  early_integration    = early_parts[, 13],
  deconvolution    = early_parts[, 14]
)   

## create ground_truth dataframe : 
gt_idx <- grep("^groundtruth", fnames_clean)

gt_names_only <- fnames_clean[gt_idx]

gt_parts <- str_split_fixed(gt_names_only, "_", 3)
gt_numbers <- str_extract(gt_parts[, 1], "\\d+$")



dataset_names <- paste0(gt_parts[, 2], gt_numbers)

dt_groundtruth <- data.table(
  full_path = files[gt_idx],
  dataset    =dataset_names
)

#####################################
#### scoring : 
########################################


library(data.table)
library(future.apply)


plan(multicore) 


all_datasets <- intersect(dt_groundtruth$dataset, late_int$dataset)




target_types <- c("basal", "classic")


score_prediction <- function(meta_row, a_real) {
  # Extract path from the metadata row
  pred_path <- meta_row$full_path 
  
  # --- Safe read attempt ---
  a_pred <- tryCatch({
    val <- read_hdf5(pred_path)$pred 
    as.matrix(val)
  }, error = function(e) {
    # Log the specific file that failed
    warning(paste("Skipping prediction file due to error:", pred_path, "-", e$message))
    return(NULL) 
  })
  
  if (is.null(a_pred)) {
    return(NULL)
  }

  a_pred_sub <- a_pred[target_types, , drop = FALSE]

  scores_full = scoring_function(A_real = a_real, A_pred = a_pred_sub)
  # scores_full = scoring_function(A_real = a_real, A_pred = a_pred)
  
  res_dt <- as.data.table(t(as.matrix(scores_full)))
  colnames(res_dt) <- names(scores_full)
  
  meta_dt <- as.data.table(meta_row)
  meta_dt[, full_path := NULL] 
  
  final_dt <- cbind(meta_dt, res_dt)
  
  return(final_dt)
}

## LI 
final_output <- future_lapply(all_datasets, function(ds) {
  
  gt_path <- dt_groundtruth[dataset == ds, full_path]
  a_real  <- read_hdf5(gt_path)$groundtruth
  a_real = as.matrix(a_real)

  a_real_sub <- a_real[target_types, , drop = FALSE]


  # B. Get all prediction paths associated with this dataset
  late_subset <- late_int[dataset == ds]
  # pred_paths <- late_int[dataset == ds, full_path]
  
  dataset_scores_list <- lapply(split(late_subset, seq_len(nrow(late_subset))), function(row) {
    score_prediction(row, a_real_sub)
  })


  return( rbindlist(dataset_scores_list) )

  
}, future.seed = TRUE)



final_results_df <- rbindlist(final_output)


datasets_string <- paste(all_datasets, collapse = "_")
file_name <- paste0("results_li_basal_classic_", datasets_string, ".csv.gz")

# file_name = "results_li_basal_classic" +  all_datasets +  ".csv.gz"
fwrite(final_results_df, file = file_name, compress = "gzip")


#########################
## ei : 
######################

final_output_ei <- future_lapply(all_datasets, function(ds) {
  
  gt_path <- dt_groundtruth[dataset == ds, full_path]
  a_real  <- read_hdf5(gt_path)$groundtruth
  a_real = as.matrix(a_real)


  a_real_sub <- a_real[target_types, , drop = FALSE]

  # B. Get all prediction paths associated with this dataset
  early_subset <- early_int[dataset == ds]
  # pred_paths <- late_int[dataset == ds, full_path]
  
  dataset_scores_list <- lapply(split(early_subset, seq_len(nrow(early_subset))), function(row) {
    score_prediction(row, a_real_sub)
  })


  return( rbindlist(dataset_scores_list) )

  
}, future.seed = TRUE)



final_results_df_ei <- rbindlist(final_output_ei)

datasets_string <- paste(all_datasets, collapse = "_")
file_name_ei <- paste0("results_ei_basal_classic_", datasets_string, ".csv.gz")

# file_name_ei = "results_li_basal_classic" +  all_datasets +  ".csv.gz"
fwrite(final_results_df_ei, file = file_name_ei, compress = "gzip")
