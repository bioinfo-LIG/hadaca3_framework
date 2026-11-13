# Supplementary Figure 3: ranking of the methods in each metric family, and size of the differences.

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

# ---- supp_rank_alternatives (alternative 2, metric families) -------------------------------------------------------
fam_cols = c(names(family_defs), "score")
fam_labels = setNames(c(names(family_defs), "Final\naggregate"), fam_cols)
fs = 6; pt = 1 / ggplot2::.pt
alt_theme = theme_bw(base_size = fs) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(size = fs, face = "bold", margin = margin(1.5, 1, 1.5, 1)),
        axis.text = element_text(size = fs), axis.text.x = element_text(angle = 45, hjust = 1, lineheight = 0.85),
        legend.position = "bottom", legend.text = element_text(size = fs), legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5))
method_scale = scale_colour_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL)
lim_diff = -0.3 # differences below are drawn at the bottom edge (MKL, MCIA, onlyRNA in some scorings)

# medians and difference to onlyDNAm, for the scorings `cols` (columns of the per-pipeline table `d`)
diff_table = function(d, cols, labels) bind_rows(lapply(cols, function(col) median_rank(d, col) %>% mutate(scoring = labels[[col]]))) %>%
  mutate(scoring = factor(scoring, levels = unname(labels[cols]))) %>%
  group_by(dataset, scoring) %>% mutate(diff = median_score - median_score[method == "OnlyMet"]) %>% ungroup()

alt1 = function(tab, title) ggplot(tab, aes(x = scoring, y = pmax(diff, lim_diff), group = method, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.4) + geom_point(size = 0.9) +
  facet_wrap(~ dataset, ncol = 4) + method_scale +
  scale_y_continuous(breaks = seq(-0.3, 0.1, 0.1), labels = function(x) sprintf("%.1f", x)) +
  labs(x = NULL, y = sprintf("median score - median score of onlyDNAm\n(cut at %.1f)", lim_diff), title = title) +
  alt_theme + guides(colour = guide_legend(nrow = 1))

summary_panel = function(tab) {
  s = tab %>% filter(method != "OnlyMet") %>% group_by(scoring, method) %>%
    summarise(med = median(diff, na.rm = TRUE), lo = min(diff, na.rm = TRUE), hi = max(diff, na.rm = TRUE), .groups = "drop") # no SDE5 restricted to 2 cell types
  ggplot(s, aes(y = method, colour = method)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_segment(aes(x = pmax(lo, lim_diff), xend = pmax(hi, lim_diff), yend = method), linewidth = 0.4) +
    geom_point(data = tab %>% filter(method != "OnlyMet", !is.na(diff)), aes(x = pmax(diff, lim_diff)), size = 0.5, alpha = 0.5) +
    geom_point(aes(x = pmax(med, lim_diff)), size = 1.6, shape = 21, fill = "white", stroke = 0.6) +
    facet_wrap(~ scoring, nrow = 1) + method_scale +
    scale_y_discrete(limits = rev(setdiff(method_order, "OnlyMet")), labels = method_label) +
    scale_x_continuous(breaks = c(-0.2, -0.1, 0), labels = function(x) sprintf("%.1f", x)) + # no label at the panel edges: they overlapped
    theme(panel.spacing.x = unit(6, "pt")) +
    labs(x = sprintf("median score - median score of onlyDNAm (cut at %.1f); point: median over the datasets, bar: range", lim_diff), y = NULL) +
    alt_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")
}

# bootstrap of the ranks: in each dataset and scoring, the pipelines of each method are resampled, medians re-ranked
boot_ranks = function(d, cols, labels, n_boot = 200) {
  set.seed(1)
  bind_rows(lapply(cols, function(col) {
    bm = d %>% filter(!is.na(.data[[col]])) %>% group_by(dataset, method) %>%
      summarise(b = list(replicate(n_boot, median(sample(.data[[col]], replace = TRUE)))), .groups = "drop") %>%
      tidyr::unnest_longer(b, indices_to = "iter") %>%
      group_by(dataset, iter) %>% mutate(rank = rank(-b, ties.method = "average")) %>%
      group_by(dataset, method) %>% summarise(lo = quantile(rank, 0.025), hi = quantile(rank, 0.975), .groups = "drop")
    bm %>% mutate(scoring = labels[[col]])
  })) %>% mutate(scoring = factor(scoring, levels = unname(labels[cols])), method = factor(method, levels = method_order))
}
alt3 = function(tab, br, title) {
  # each method shifted by a small fixed offset, so that its interval sits on its point
  dd = tab %>% left_join(br, by = c("dataset", "method", "scoring")) %>%
    mutate(x = as.integer(scoring) + (as.integer(method) - (length(method_order) + 1) / 2) * 0.05)
  ggplot(dd, aes(x = x, y = rank, group = method, colour = method)) +
    geom_line(linewidth = 0.35, alpha = 0.7) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6) +
    geom_point(size = 0.8) +
    scale_x_continuous(breaks = seq_along(levels(tab$scoring)), labels = levels(tab$scoring)) +
    scale_y_reverse(breaks = seq_along(method_order)) +
    facet_wrap(~ dataset, ncol = 4) + method_scale +
    labs(x = NULL, y = "Median rank (1 = best), 95 % bootstrap interval", title = title) +
    alt_theme + guides(colour = guide_legend(nrow = 1))
}
bump_small = function(tab, title) bump_plot(tab, title) + alt_theme + guides(colour = guide_legend(nrow = 1)) +
  geom_line(linewidth = 0.4) + geom_point(size = 0.9)

tab = diff_table(pipelines_noVIVO, fam_cols, fam_labels)
p2 = wrap_plots(bump_small(tab, "Metric families and final aggregate score") + labs(tag = "A"), summary_panel(tab) + labs(tag = "B"), ncol = 1, heights = c(1.6, 1)) &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(fig_dir, "SuppFig3.pdf"), plot = p2, width = 18, height = 18, units = "cm", device = cairo_pdf)
