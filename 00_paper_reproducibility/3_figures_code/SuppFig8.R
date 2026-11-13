# Supplementary Figure 8: best pipeline of each method compared with its uni-omic counterparts.

data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_03_figure5.R")

# ---- rows ----------------------------------------------------------------------------------------------------------
int_methods = c(early_methods, late_methods) # rows of panels A and B (first method on top)
y_int = function(m) length(int_methods) + 1 - match(as.character(m), int_methods) # numeric row position
regime_sep = geom_hline(yintercept = length(late_methods) + 0.5, colour = "grey40", linewidth = 0.4) # early above, late below
y_scale = scale_y_continuous(breaks = y_int(int_methods), labels = method_label(int_methods), expand = expansion(add = 0.5))

# ---- supp_top_vs_uniomic_data --------------------------------------------------------------------------------------
shown_methods = int_methods
top_pipelines = bind_rows(ei_ds %>% mutate(regime = "Early"),
                          li_ds %>% filter(!method %in% uniomic_methods) %>% mutate(regime = "Late")) %>%
  group_by(regime, dataset, method) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(regime, dataset, method, score, key_met, key_rna) %>%
  left_join(select(onlymet, key_met, score_matched_dnam = score_onlymet), by = "key_met") %>%
  left_join(select(onlyrna, key_rna, score_matched_rna = score_onlyrna), by = "key_rna") %>%
  left_join(best_dna, by = "dataset") %>%
  left_join(best_rna, by = "dataset") %>%
  mutate(method = factor(method, levels = method_order), dataset_lab = lab_dataset(dataset),
         gain_vs_matched_dnam = score - score_matched_dnam,
         gain_vs_best_dnam = score - best_onlyDNAm, gain_vs_best_rna = score - best_onlyRNA) %>%
  arrange(dataset, method)

lab_int = "integrated top pipeline"
lab_dna = "onlyDNAm, same DNAm side"
lab_rna = "onlyRNA, same RNA side"

best_lines = bind_rows(
  best_dna %>% transmute(dataset, which = "best onlyDNAm of the dataset", x = best_onlyDNAm),
  best_rna %>% transmute(dataset, which = "best onlyRNA of the dataset", x = best_onlyRNA)
) %>% mutate(dataset_lab = lab_dataset(dataset))

# ---- supp_top_vs_uniomic -------------------------------------------------------------------------------------------
xmin_plot = 0.45
xmax_plot = 1.0
band_start = 0.935

# for the same `size`, a square (shape 22) looks smaller than a triangle (shape 24): scaled to look equally big
size_tri = 2.1
size_sq = 1.4 * size_tri
stroke_small = 0.7

dumb = top_pipelines %>%
  mutate(rna_off = score_matched_rna < xmin_plot, rna_plot = pmax(score_matched_rna, xmin_plot))

p_dumbbell = ggplot() +
  annotate("rect", xmin = band_start, xmax = xmax_plot, ymin = -Inf, ymax = Inf, fill = "grey95") +
  geom_vline(data = best_lines, aes(xintercept = x, linetype = which), colour = "grey35", linewidth = 0.5) +
  geom_segment(data = dumb, aes(x = score_matched_dnam, xend = score, y = method, yend = method),
               colour = "grey70", linewidth = 0.8) +
  geom_point(data = filter(dumb, !rna_off), aes(x = rna_plot, y = method, shape = lab_rna),
             fill = "grey85", colour = "grey60", size = size_sq, stroke = stroke_small) +
  geom_point(data = filter(dumb, rna_off), aes(x = rna_plot, y = method, shape = lab_rna),
             fill = "white", colour = "firebrick", size = size_sq, stroke = stroke_small) +
  geom_text(data = filter(dumb, rna_off), aes(x = rna_plot, y = method, label = sprintf("%.2f", score_matched_rna)),
            colour = "firebrick", hjust = -0.35, size = 6 / .pt) +
  geom_point(data = dumb, aes(x = score_matched_dnam, y = method, shape = lab_dna),
             fill = "white", colour = "grey20", size = size_tri, stroke = stroke_small) +
  geom_point(data = dumb, aes(x = score, y = method, fill = method, shape = lab_int),
             colour = "black", size = 3.6, stroke = 0.6) +
  geom_text(data = dumb, aes(x = xmax_plot - 0.006, y = method, label = sprintf("%.2f", score)),
            hjust = 1, size = 6.5 / .pt, fontface = "bold") +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_shape_manual(values = setNames(c(21, 24, 22), c(lab_int, lab_dna, lab_rna)),
                     breaks = c(lab_int, lab_dna, lab_rna), name = NULL) +
  scale_linetype_manual(values = c("best onlyDNAm of the dataset" = "solid", "best onlyRNA of the dataset" = "dashed"), name = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(
    fill = c("grey80", "white", "grey85"), colour = c("black", "grey20", "grey60"),
    size = c(3.6, size_tri, size_sq), stroke = c(0.6, stroke_small, stroke_small)))) +
  scale_x_continuous(breaks = seq(0.5, 0.9, 0.1)) +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot)) +
  scale_y_discrete(limits = rev(shown_methods), labels = method_label) +
  facet_wrap(~ dataset_lab, ncol = 3) +
  labs(x = "Aggregated score (number on the right: score of the integrated top pipeline)", y = element_blank()) +
  theme_bw(base_size = 7) +
  theme(
    strip.text = element_text(size = 7), axis.text = element_text(size = 6.5), axis.title.x = element_text(size = 7),
    panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey92"),
    legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 7),
    legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0), legend.key.height = unit(9, "pt")
  )

p_dumbbell = p_dumbbell + labs(title = "Score of the best pipeline vs. its uni-omic counterparts") + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8.5))
# A4 page: 18 cm wide
ggsave(file.path(fig_dir, "SuppFig8.pdf"), plot = p_dumbbell, width = 18, height = 20, units = "cm")
