# Figure 5 of the article: top-performing pipelines.

data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_03_figure5.R")

# ---- rows ----------------------------------------------------------------------------------------------------------
int_methods = c(early_methods, late_methods) # rows of panels A and B (first method on top)
y_int = function(m) length(int_methods) + 1 - match(as.character(m), int_methods) # numeric row position
regime_sep = geom_hline(yintercept = length(late_methods) + 0.5, colour = "grey40", linewidth = 0.4) # early above, late below
y_scale = scale_y_continuous(breaks = y_int(int_methods), labels = method_label(int_methods), expand = expansion(add = 0.5))

# ---- panel_A_data --------------------------------------------------------------------------------------------------
peak = bind_rows(ei_ds, li_ds %>% filter(method %in% late_methods)) %>%
  group_by(dataset, method) %>%
  summarise(n_pipelines = n(), top_score = max(score), .groups = "drop") %>%
  left_join(best_dna, by = "dataset") %>%
  mutate(diff = top_score - best_onlyDNAm) %>%
  group_by(dataset) %>% mutate(best_of_dataset = diff == max(diff)) %>% ungroup()

# ---- gain_summary --------------------------------------------------------------------------------------------------
top_frac = 0.1 # the "top" of a (dataset, method) group = its top 10 % pipelines
gain_summary = bind_rows(ei_ds, li_ds %>% filter(method %in% late_methods)) %>%
  left_join(select(onlymet, key_met, score_onlymet), by = "key_met") %>%
  filter(!is.na(score_onlymet)) %>%
  mutate(delta_met = score - score_onlymet) %>%
  group_by(dataset, method) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  summarise(n_pipelines = n(), n_top = ceiling(top_frac * n()),
            delta_all = median(delta_met), delta_top = median(delta_met[seq_len(first(n_top))]), .groups = "drop")

# ---- panel_A_rank_data ---------------------------------------------------------------------------------------------
top_scores = bind_rows(ei_ds, li_ds) %>%
  group_by(dataset, method) %>%
  summarise(top_score = max(score), .groups = "drop") %>%
  group_by(dataset) %>% mutate(rank = rank(-top_score, ties.method = "min")) %>% ungroup() %>%
  mutate(method = factor(method, levels = method_order), dataset = factor(as.character(dataset), levels = dataset_levels))

# ---- figure5 -------------------------------------------------------------------------------------------------------
pt = 1 / ggplot2::.pt
tag5 = theme(plot.tag = element_text(face = "bold", size = 14))
# method colours of all the figures (pale shades): points with a dark outline, so that the palest shades stay visible

# ---- A: rank of the top score ----
p_A = ggplot(top_scores, aes(x = as.integer(dataset), y = rank, colour = method, group = method)) +
  geom_vline(xintercept = length(dataset_levels) - 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_line(data = filter(top_scores, dataset != "VIVO"), linewidth = 0.55) +
  # as figure 3D: dotted link from the last in silico dataset to VIVO, VIVO points as triangles
  geom_line(data = filter(top_scores, as.integer(dataset) >= length(dataset_levels) - 1), linewidth = 0.45, linetype = "dotted") +
  geom_point(aes(fill = method, shape = dataset == "VIVO"), colour = "grey30", size = 1.4, stroke = 0.3) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
  scale_colour_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL, guide = "none") +
  scale_fill_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL) +
  scale_x_continuous(breaks = seq_along(dataset_levels), labels = ifelse(dataset_levels == "VIVO", "VIVO*", dataset_levels), expand = expansion(add = 0.4)) +
  scale_y_reverse(breaks = 1:10) +
  labs(title = "Rank of the top score", x = NULL, y = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1), axis.text.y = element_text(size = 6),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5), legend.position = "bottom", legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"), legend.key.width = unit(0.4, "cm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) +
  guides(fill = guide_legend(ncol = 3, override.aes = list(shape = 22, size = 3, colour = NA, stroke = 0))) # squares, as figure 3B; under A: more width for the heatmap B (colour and fill legends merged)

# ---- B: peak score minus the best onlyDNAm (heatmap) ----
lim5 = 0.05 # colour limit (squished beyond: the values are written in the cells)
dB = peak %>% mutate(dataset = factor(as.character(dataset), levels = dataset_levels),
                     block = factor(ifelse(dataset == "VIVO", "VIVO*", "benchmark"), levels = c("benchmark", "VIVO*")),
                     y = y_int(method), lab = ifelse(abs(diff) < 0.0005, "0.000", sprintf("%+.3f", diff)))
p_B = ggplot(dB, aes(x = dataset, y = y, fill = diff)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = filter(dB, best_of_dataset), fill = NA, colour = "black", linewidth = 0.6, width = 0.94, height = 0.92) +
  geom_text(aes(label = lab, colour = abs(diff) > 0.75 * lim5), size = 6 * pt) +
  regime_sep +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                       limits = c(-lim5, lim5), oob = scales::squish, name = "Score difference",
                       guide = guide_colourbar(barwidth = unit(3.5, "cm"), barheight = unit(0.22, "cm"), title.position = "top", title.hjust = 0.5)) +
  y_scale +
  facet_grid(. ~ block, scales = "free_x", space = "free_x") +
  labs(title = "Peak: best integrated - best onlyDNAm", x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), panel.spacing.x = unit(4, "pt"), axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.5), strip.text = element_text(size = 6.5, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5), legend.position = "bottom",
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0))

# ---- C: score gain from adding the RNA, one panel per dataset ----
xlim_C = c(-0.35, 0.3) # values below: drawn at the edge, with their value (Kernel, all the pipelines)
lab_all = "all pipelines"; lab_top = "top 10 % pipelines"
dC = bind_rows(gain_summary %>% transmute(dataset, method, which = lab_all, delta = delta_all),
               gain_summary %>% transmute(dataset, method, which = lab_top, delta = delta_top)) %>%
  mutate(which = factor(which, levels = c(lab_all, lab_top)), off = delta < xlim_C[1], x = pmax(delta, xlim_C[1]),
         method = factor(method, levels = int_methods),
         dataset = factor(ifelse(dataset == "VIVO", "VIVO*", as.character(dataset)), levels = c(setdiff(dataset_levels, "VIVO"), "VIVO*")))
segC = gain_summary %>% transmute(method = factor(method, levels = int_methods), x0 = pmax(delta_all, xlim_C[1]), x1 = pmax(delta_top, xlim_C[1]),
                                  dataset = factor(ifelse(dataset == "VIVO", "VIVO*", as.character(dataset)), levels = levels(dC$dataset)))
p_C = ggplot(dC, aes(x = x, y = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = length(late_methods) + 0.5, colour = "grey60", linewidth = 0.3) + # between uniPort and limean (reversed y: late methods at the bottom)
  geom_segment(data = segC, aes(x = x0, xend = x1, y = method, yend = method), colour = "grey75", linewidth = 0.6) +
  geom_point(data = filter(dC, which == lab_all), aes(shape = which), fill = "white", colour = "grey40", size = 1.4, stroke = 0.5) +
  geom_point(data = filter(dC, which == lab_top), aes(shape = which, fill = method), colour = "grey25", size = 1.5, stroke = 0.3) +
  geom_text(data = filter(dC, off), aes(label = sprintf("%.2f", delta)), hjust = 0, vjust = -0.7, size = 6 * pt, colour = "grey30") + # above the point
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_shape_manual(values = setNames(c(21, 21), c(lab_all, lab_top)), name = NULL,
                     guide = guide_legend(override.aes = list(fill = c("white", "grey60"), colour = c("grey40", "grey25")))) +
  scale_y_discrete(limits = rev(int_methods), labels = method_label) +
  scale_x_continuous(breaks = c(-0.3, -0.15, 0, 0.15, 0.3)) +
  coord_cartesian(xlim = xlim_C, clip = "off") +
  facet_wrap(~ dataset, ncol = 3) +
  labs(title = "RNA gain, same DNAm side", x = "integrated pipeline - onlyDNAm pipeline with the same DNAm side (median)", y = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey93"), strip.text = element_text(size = 6.5, face = "bold"),
        axis.text = element_text(size = 6), axis.title.x = element_text(size = 6.5), plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5),
        legend.position = "bottom", legend.text = element_text(size = 6.5), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))

figure5 = wrap_plots(wrap_plots(p_A + labs(tag = "A") + tag5, p_B + labs(tag = "B") + tag5, nrow = 1, widths = c(0.62, 1)),
                     p_C + labs(tag = "C") + tag5, ncol = 1, heights = c(0.8, 1.25))

# ---- figure5_split -------------------------------------------------------------------------------------------------
p_A_side = p_A + labs(tag = "A", y = "Rank (1 = best)") + tag5 + theme(legend.position = "none")
# legend of the methods: squares, on one line under the whole figure (as figure 3B)
dev_legend = cairo_pdf(tempfile(fileext = ".pdf")) # the legend is measured with the device of the final figure (same result in the Rmd and in the scripts)
legend_A = ggplotGrob(p_A + guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 22, size = 3, colour = NA, stroke = 0))) +
                        theme(legend.position = "bottom"))
legend_A = legend_A$grobs[[which(legend_A$layout$name == "guide-box-bottom")]]
dev_legend = dev.off()
figure5_main = wrap_plots(wrap_plots(p_A_side, p_C + labs(tag = "B") + tag5, nrow = 1, widths = unit(c(6.7, 1), c("cm", "null"))), # panel A as wide as figure 3D (6.7 cm)
                          wrap_elements(full = legend_A), ncol = 1, heights = c(1, 0.06))
ggsave(file.path(fig_dir, "Fig5.pdf"), plot = figure5_main, width = 18, height = 12, units = "cm", device = cairo_pdf)
