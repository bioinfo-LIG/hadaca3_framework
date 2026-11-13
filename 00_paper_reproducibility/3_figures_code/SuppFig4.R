# Supplementary Figure 4: radar plots of the metrics of the best pipeline of each method.

data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_01_figure3.R")

# ---- supp_data -----------------------------------------------------------------------------------------------------
ei_ds = read_res(file.path(data_dir, "ei_ds.rds"))
li_ds = read_res(file.path(data_dir, "li_ds.rds"))
ds_noVIVO = setdiff(dataset_levels, "VIVO")

# the 4 metric families of the aggregated score: geometric mean of their normalised metrics. A family is computed on
# the metrics available (JSD is missing in VITR and SDN4: their error family uses RMSE, MAE and Aitchison)
family_defs = list(
  "Cohort\ncorrelation" = c("pearson_tot_norm", "spearman_tot_norm", "pearson_col_norm", "spearman_col_norm"),
  "Cell-type\ncorrelation" = c("pearson_row_norm", "spearman_row_norm"),
  "Error\nmetrics" = c("rmse_norm", "mae_norm", "aitchison_norm", "jsd_norm"),
  "Geometric\nmetrics" = c("aid_norm", "sdid_norm"))
geo_mean = function(m) exp(rowMeans(log(m), na.rm = TRUE))
pipelines_noVIVO = bind_rows(ei_ds, li_ds) %>% filter(dataset %in% ds_noVIVO) %>% mutate(dataset = factor(as.character(dataset), levels = ds_noVIVO))
for (f in names(family_defs)) pipelines_noVIVO[[f]] = geo_mean(as.matrix(pipelines_noVIVO[family_defs[[f]]]))
# check: the geometric mean of the 4 families is the stored aggregated score (exact, except in VITR and SDN4, where
# JSD is missing: the stored score then drops the weight of JSD instead of averaging the error family on 3 metrics)
check = geo_mean(as.matrix(pipelines_noVIVO[names(family_defs)]))

# median score of each method in each dataset, and rank of the methods in the dataset (1 = best)
median_rank = function(d, col) d %>% group_by(dataset, method) %>%
  summarise(median_score = median(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
  group_by(dataset) %>% mutate(rank = rank(-median_score, ties.method = "average", na.last = "keep")) %>% ungroup() %>% # no score (SDE5 restricted to 2 cell types): no rank
  mutate(method = factor(method, levels = method_order))

# bump chart: rank of each method (lines) for each scoring (x), one panel per dataset
bump_plot = function(d, title) {
  ggplot(d, aes(x = scoring, y = rank, group = method, colour = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    scale_y_reverse(breaks = seq_along(method_order)) +
    scale_colour_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL) +
    facet_wrap(~ dataset, ncol = 4) +
    labs(x = NULL, y = "Median rank (1 = best)", title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, lineheight = 0.9), strip.text = element_text(size = 8),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    guides(colour = guide_legend(nrow = 1))
}

# ---- supp_figure5 --------------------------------------------------------------------------------------------------
radar_metrics = c(aid_norm = "aid", sdid_norm = "sdid",
                  pearson_tot_norm = "pr_t", spearman_tot_norm = "sp_t", pearson_col_norm = "pr_c", spearman_col_norm = "sp_c",
                  pearson_row_norm = "pr_r", spearman_row_norm = "sp_r",
                  rmse_norm = "rmse", mae_norm = "mae", aitchison_norm = "aitch.", jsd_norm = "jsd")
radar_family = rep(c("Geometric", "Cohort correlation", "Cell-type correlation", "Error"), c(2, 4, 2, 4))
family_text_cols = c("Geometric" = "#E41A1C", "Cohort correlation" = "#C266C6", "Cell-type correlation" = "#F28E2B", "Error" = "#3F4AA8")
dataset_cols = setNames(brewer.pal(8, "Dark2"), ds_noVIVO)

# best pipeline of each method: best mean score over the 9 datasets (pipelines evaluated on the 9 datasets), as the
# top combinations of the supplementary tables; its metrics are shown on the 8 datasets without VIVO
id_cols = union(all_functions_ei, all_functions_li)
make_pipeline_id = function(d) do.call(paste, c(lapply(d[intersect(id_cols, names(d))], as.character), sep = "|"))
best_ids = bind_rows(ei_ds, li_ds) %>%
  mutate(id = make_pipeline_id(.)) %>%
  group_by(method, id) %>% filter(n_distinct(dataset) == length(dataset_levels)) %>%
  summarise(mean_score = mean(score), .groups = "drop") %>%
  group_by(method) %>% slice_max(mean_score, n = 1, with_ties = FALSE) %>% ungroup()
best_pipelines = pipelines_noVIVO %>% mutate(id = make_pipeline_id(.)) %>%
  inner_join(select(best_ids, method, id, mean_score), by = c("method", "id"))

n_m = length(radar_metrics)
angle = pi / 2 - 2 * pi * (seq_len(n_m) - 1) / n_m          # first metric on top, clockwise
grid = bind_rows(lapply(c(0.2, 0.4, 0.6, 0.8, 1), function(r) data.frame(r = r, k = c(seq_len(n_m), 1)))) %>%
  mutate(x = r * cos(angle[k]), y = r * sin(angle[k]))
spokes = data.frame(x = 0, y = 0, xend = cos(angle), yend = sin(angle))
lab = data.frame(label = radar_metrics, family = radar_family, x = 1.28 * cos(angle), y = 1.22 * sin(angle))

radar_plot = function(m, group = NULL) {
  d = best_pipelines %>% filter(method == m) %>% select(dataset, all_of(names(radar_metrics)))
  poly = bind_rows(lapply(seq_len(nrow(d)), function(i) {
    v = unlist(d[i, names(radar_metrics)])
    data.frame(dataset = d$dataset[i], pos = seq_len(n_m + 1), k = c(seq_len(n_m), 1), r = c(v, v[1]))
  })) %>% mutate(x = r * cos(angle[k]), y = r * sin(angle[k]))
  # polygon closed over a missing metric (JSD in VITR and SDN4): dashed segment between the two metrics around it
  bridge = poly %>% filter(!is.na(r)) %>% group_by(dataset) %>%
    mutate(xend = lead(x), yend = lead(y), skip = lead(pos) - pos > 1) %>% ungroup() %>% filter(skip %in% TRUE)
  ggplot() +
    geom_path(data = grid, aes(x, y, group = r), colour = "grey80", linewidth = 0.25) +
    geom_segment(data = spokes, aes(x = x, y = y, xend = xend, yend = yend), colour = "grey80", linewidth = 0.25) +
    geom_path(data = poly, aes(x, y, colour = dataset, group = dataset), linewidth = 0.45, na.rm = TRUE) +
    geom_segment(data = bridge, aes(x = x, y = y, xend = xend, yend = yend, colour = dataset), linewidth = 0.45, linetype = "22") +
    geom_point(data = poly, aes(x, y, colour = dataset), size = 0.4, na.rm = TRUE) +
    geom_text(data = lab, aes(x, y, label = label), colour = family_text_cols[lab$family], size = 6 / .pt, fontface = "bold") +
    scale_colour_manual(values = dataset_cols, name = NULL, drop = FALSE) +
    coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.35, 1.35), clip = "off") +
    labs(title = method_label(m), subtitle = group) +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 6.5, colour = "grey30"),
          legend.position = "bottom", legend.text = element_text(size = 7)) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8, linewidth = 0.8)))
}

# legend of the metric families, as text coloured by family
fam_text = data.frame(
  family = rep(c("Geometric", "Cohort correlation", "Error", "Cell-type correlation"), c(3, 5, 5, 3)),
  text = c("Geometric family", "aid: angular inequality/disproportionality (ID)", "sdid: sine-diagonal ID",
           "Cohort correlation - intra-sample variations", "pr_t: total Pearson correlation", "sp_t: total Spearman correlation",
           "pr_c: column Pearson correlation", "sp_c: column Spearman correlation",
           "Error family", "rmse: root mean square error", "mae: mean absolute error", "aitch.: Aitchison distance", "jsd: Jensen-Shannon divergence",
           "Cell-type correlation - inter-sample variations", "pr_r: cell-type Pearson correlation", "sp_r: cell-type Spearman correlation"),
  head = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)) %>%
  mutate(y = -seq_along(text) - cumsum(head) * 0.4) # one column (A4 page: the two-column layout overlapped)
p_fam = ggplot(fam_text, aes(x = 0, y = y, label = text, colour = family)) +
  geom_text(aes(fontface = ifelse(head, "bold", "plain")), hjust = 0, size = 6 / .pt) +
  scale_colour_manual(values = family_text_cols, guide = "none") +
  coord_cartesian(xlim = c(0, 1), clip = "off") + theme_void()

# A4 page (18 cm wide): 4 radars per row; the type of each method is given under its name
plots = c(lapply(early_methods, radar_plot, group = "early integration"),
          lapply(late_methods, radar_plot, group = "late integration"),
          lapply(uniomic_methods, radar_plot, group = "uni-omic"),
          list(p_fam))
design = "ABCD\nEFGH\nIJKK"
supp_figure5 = wrap_plots(plots, design = design, guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "SuppFig4.pdf"), plot = supp_figure5, width = 18, height = 19, units = "cm")
