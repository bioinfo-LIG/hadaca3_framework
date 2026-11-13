

# ---- re_loading_pckgs ----------------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
theme_set(theme_light())
library(RColorBrewer)
library(wesanderson)
library(patchwork)

# ---- settings ------------------------------------------------------------------------------------------------------
out_dir = tempdir() # intermediate files of the analysis (not kept)

early_methods = c("concatnoscale", "concatscale", "Kernel", "omicade4bulk", "uniPort")
late_methods = c("limean", "limeanRMSE", "TunedJ")
uniomic_methods = c("OnlyMet", "OnlyRna")
method_order = c(early_methods, late_methods, uniomic_methods)

modules_early = list("PP RNA" = "preprocessing_RNA", "PP DNAm" = "preprocessing_MET",
                     "FS RNA" = "feature_selection_RNA", "FS DNAm" = "feature_selection_MET",
                     "Deconv" = "deconvolution")
modules_late = list("PP RNA" = "preprocessing_RNA", "PP DNAm" = "preprocessing_MET",
                    "FS RNA" = "feature_selection_RNA", "FS DNAm" = "feature_selection_MET",
                    "Deconv RNA" = "deconvolution_rna", "Deconv DNAm" = "deconvolution_met")

all_functions_ei = c('preprocessing_mixRNA', 'feature_selection_mixRNA', 'preprocessing_RNA', 'feature_selection_RNA',
                     'preprocessing_scRNA', 'feature_selection_scRNA', 'preprocessing_mixMET', 'feature_selection_mixMET',
                     'preprocessing_MET', 'feature_selection_MET', 'early_integration', 'deconvolution')
all_functions_li = c('preprocessing_mixRNA', 'feature_selection_mixRNA', 'preprocessing_RNA', 'feature_selection_RNA',
                     'preprocessing_scRNA', 'feature_selection_scRNA', 'deconvolution_rna', 'preprocessing_mixMET',
                     'feature_selection_mixMET', 'preprocessing_MET', 'feature_selection_MET', 'deconvolution_met', 'late_integration')

# same levels as the chunks ei_levels_def and li_levels_def of the paper code
lev_ei = list(
  preprocessing_RNA = c("ppID", "LogNorm", "SumNorm"),
  feature_selection_RNA = c("fsID", "SCcluster", "scpseudobulk", "Toastbulknbfs", "Toastvst"),
  preprocessing_scRNA = c("ppID", "sccluster", "scConcat", "LogNorm", "SumNorm", "scbinarypseudobulklog"),
  feature_selection_scRNA = c("fsID", "SCcluster", "scpseudobulk", "Toastbulknbfs", "Toastvst"),
  preprocessing_MET = c("ppID", "LogNorm", "SumNorm"),
  feature_selection_MET = c("fsID", "maxdiscriminant", "mostmethylated", "SPLSDA", "Toastpercent"),
  early_integration = early_methods,
  deconvolution = c("nnls", "lm", "RLR", "RLRnnls", "wRLR", "epic")
)
lev_li = list(
  preprocessing_RNA = c("ppID", "LogNorm", "SumNorm", "nopp"),
  feature_selection_RNA = c("fsID", "SCcluster", "scpseudobulk", "Toastbulknbfs", "Toastvst", "nofs"),
  preprocessing_scRNA = c("ppID", "LogNorm", "SumNorm", "scbinarypseudobulklog", "sccluster", "scConcat", "nopp"),
  feature_selection_scRNA = c("fsID", "SCcluster", "scpseudobulk", "Toastbulknbfs", "Toastvst", "nofs"),
  preprocessing_MET = c("ppID", "LogNorm", "SumNorm", "nopp"),
  feature_selection_MET = c("fsID", "maxdiscriminant", "mostmethylated", "SPLSDA", "Toastpercent", "nofs"),
  deconvolution_rna = c("nnls", "lm", "RLR", "RLRnnls", "wRLR", "epic", "node"),
  deconvolution_met = c("nnls", "lm", "RLR", "RLRnnls", "wRLR", "epic", "node"),
  late_integration = c("limean", "TunedJ", "limeanRMSE", "OnlyMet", "OnlyRna")
)
# datasets: VIVO last (its aggregated score uses a different formula)
dataset_levels = c("VITR", "SBN5", "SDN5", "SDN4", "SDN6", "SDE5", "SDEL", "SDC5", "VIVO")
dataset_labels = c(setdiff(dataset_levels, "VIVO"), "VIVO (different score)")
lab_dataset = function(d) factor(ifelse(d == "VIVO", "VIVO (different score)", as.character(d)), levels = dataset_labels)

# ---- color_choices -------------------------------------------------------------------------------------------------
# ---- modules: one colour per module family ----
# sand / ochre / copper for the modules, two olives for the integration: none of the colours of the methods (violet, turquoise,
# grey) or of the heatmaps (red, blue)
module_pal = c(preprocessing = "#EEEAE4", feature_selection = "#DAD3C9", deconvolution = "#BDB2A3", # neutral greiges
               early_integration = "#D6EAF4", late_integration = "#FBE6C2") # pale tints of the early / late method families
module_family = c("PP RNA" = "preprocessing", "PP DNAm" = "preprocessing",
                  "FS RNA" = "feature_selection", "FS DNAm" = "feature_selection",
                  "Deconv" = "deconvolution", "Deconv RNA" = "deconvolution", "Deconv DNAm" = "deconvolution")

# ---- module levels ----
all_levels = unique(as.character(unlist(c(lev_ei, lev_li))))
id_levels = all_levels[grepl("ID", all_levels) | grepl("nopp", all_levels) | grepl("nofs", all_levels) |
                         grepl("node", all_levels) | grepl("OnlyRna", all_levels) | grepl("OnlyMet", all_levels)]
other_levels = all_levels[!(grepl("ID", all_levels) | grepl("nopp", all_levels) |
                              grepl("nofs", all_levels) | grepl("node", all_levels))]
palette_colors = colorRampPalette(c(brewer.pal(7, "Set2"), brewer.pal(7, "Dark2"), brewer.pal(8, "Set3")))(length(other_levels))
level_colors = c(setNames(rep("grey60", length(id_levels)), id_levels),
                 setNames(palette_colors, other_levels))
level_levels = unique(c(lev_ei$preprocessing_RNA, lev_ei$feature_selection_RNA,
                        lev_ei$feature_selection_MET, lev_ei$deconvolution))

# ---- integration methods: violets (early) / turquoises (late), greys for the uni-omic ----
method_colors = c(setNames(c("#BFE3F0", "#7CC4E0", "#3A9BCB", "#1F6DA8", "#103F73"), early_methods), # early integration: blues (sky to navy)
                  setNames(c("#F8D27A", "#EE9A3A", "#B9581C"), late_methods),     # late integration: ambers (yellow to rust)
                  setNames(c("grey35", "grey65"), uniomic_methods))
# names of the methods in the figures (the framework and the result tables keep the original names)
method_names = c(concatnoscale = "Concat", concatscale = "ConcatStd", Kernel = "MKL", omicade4bulk = "MCIA", uniPort = "uniPort",
                 limean = "Mean", limeanRMSE = "wMean", TunedJ = "PriorMean", OnlyMet = "onlyDNAm", OnlyRna = "onlyRNA")
method_label = function(x) { x = as.character(x); ifelse(x %in% names(method_names), method_names[x], x) }

# ---- diverging scale of the heatmaps ----
div_cols = c(low = "#762A83", mid = "#F4F3F0", high = "#1B7837") # purple (negative) / off-white (0) / green (positive)
fill_scale = function(lim) {
  scale_fill_gradient2(
    low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
    limits = c(-lim, lim), oob = scales::squish, na.value = "grey90",
    name = "Diff median aggr score",
    guide = guide_colourbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"),
                            title.position = "top", title.hjust = 0.5)
  )
}
tag_theme = theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.01, 0.99))

# ---- data ----------------------------------------------------------------------------------------------------------
# names of the figures for two blocks of the framework: the "Scale" pre-processing divides each sample by its sum (no
# standardisation): "SumNorm"; the "RLRpoisson" deconvolution is a robust regression with prior weights 1 / mean of the
# feature in the reference (not a Poisson regression): "wRLR"
new_names = list(preprocessing = c(Scale = "SumNorm"), deconvolution = c(RLRpoisson = "wRLR"))
rename_pp = function(d) {
  for (type in names(new_names)) for (cl in grep(paste0("^", type), names(d), value = TRUE)) {
    v = d[[cl]]
    for (o in names(new_names[[type]])) {
      if (is.factor(v)) levels(v)[levels(v) == o] = new_names[[type]][[o]] else v[v == o] = new_names[[type]][[o]]
    }
    d[[cl]] = v
  }
  d
}
read_res = function(f) rename_pp(readRDS(f))
res_with = list(ei = read_res(file.path(data_dir, "results_ei_withVIVO.rds")), # one row per pipeline (mean over the 9 datasets)
                li = read_res(file.path(data_dir, "results_li_withVIVO.rds")))
ei_ds = read_res(file.path(data_dir, "ei_ds.rds")) # one row per pipeline and per dataset
li_ds = read_res(file.path(data_dir, "li_ds.rds"))
all_ds = bind_rows(ei_ds %>% select(dataset, method, score, all_of(c("preprocessing_mixMET", "feature_selection_mixMET", "preprocessing_MET", "feature_selection_MET"))),
                   li_ds %>% select(dataset, method, score, all_of(c("preprocessing_mixMET", "feature_selection_mixMET", "preprocessing_MET", "feature_selection_MET"))))
rm(ei_ds, li_ds)
