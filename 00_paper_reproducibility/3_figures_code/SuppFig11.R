# Supplementary Figure 11: generalist pipeline of each strategy compared with the generalist onlyDNAm pipeline.

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

# ---- figure5_heatmap_v2 --------------------------------------------------------------------------------------------
lim_v2 = 0.04
heat_v2_cols = c(div_cols[["low"]], "#C2A5CF", div_cols[["mid"]], "#A6DBA0", div_cols[["high"]]) # off-white centre: the cells near 0 stay visible
fmt_v2 = function(v) ifelse(abs(v) < 0.0005, "0", paste0(ifelse(v > 0, "+", "-"), sub("^0", "", sprintf("%.3f", abs(v)))))
ds_bench = setdiff(dataset_levels, "VIVO")
col_x = c(setNames(seq_along(ds_bench), ds_bench), median = length(ds_bench) + 1.55, VIVO = length(ds_bench) + 2.75)
row_y = function(m) { i = match(as.character(m), int_methods); -(i + ifelse(i > length(early_methods), 0.45, 0)) } # gap between early and late
cells = bind_rows(
  dB %>% transmute(method = as.character(method), col = as.character(dataset), v = diff, win = best_of_dataset & dataset != "VIVO"), # VIVO: faded, not framed
  dB %>% filter(dataset != "VIVO") %>% group_by(method) %>% summarise(v = median(diff), .groups = "drop") %>%
    transmute(method = as.character(method), col = "median", v, win = FALSE)) %>%
  mutate(x = col_x[col], y = row_y(method), faded = col == "VIVO", bold = win | col == "median",
         ink = ifelse(faded, "grey35", ifelse(abs(v) >= 0.7 * lim_v2, "white", "grey15")))
half = 0.46 # half side of a cell (gap 0.08)
rect_poly = function(d, h = half) bind_rows(lapply(seq_len(nrow(d)), function(k)
  data.frame(id = k, px = d$x[k] + c(-h, h, h, -h), py = d$y[k] + c(-h, -h, h, h), v = d$v[k], faded = d$faded[k])))
poly = rect_poly(cells)
poly_win = rect_poly(filter(cells, win))
lab_rows = data.frame(method = int_methods) %>% mutate(y = row_y(method))
heads = data.frame(x = col_x, label = c(ds_bench, "median", "VIVO*"), face = c(rep("plain", length(ds_bench)), "bold", "plain"),
                   colour = c(rep("grey15", length(ds_bench) + 1), "grey45"))
groups = data.frame(y = c(row_y(early_methods[1]), row_y(late_methods[1])) + 0.74, label = c("Early integration", "Late integration"),
                    colour = unname(method_colors[c(early_methods[4], late_methods[3])]))
x_name = 0.35 # right edge of the method names
p_B2 = ggplot() +
  ggforce::geom_shape(data = poly, aes(x = px, y = py, group = id, fill = v, alpha = faded), radius = unit(1.2, "pt")) +
  ggforce::geom_shape(data = poly_win, aes(x = px, y = py, group = id), fill = NA, colour = "grey10", linewidth = 0.5, radius = unit(1.2, "pt")) +
  geom_text(data = cells, aes(x = x, y = y, label = fmt_v2(v), colour = ink, fontface = ifelse(bold, "bold", "plain")), size = 6 * pt) +
  geom_text(data = heads, aes(x = x, y = 0.12, label = label, fontface = face, colour = colour), size = 6.3 * pt, vjust = 0) +
  geom_text(data = groups, aes(x = -2.05, y = y, label = label, colour = colour), size = 6.5 * pt, fontface = "bold", hjust = 0) +
  geom_text(data = lab_rows, aes(x = x_name, y = y, label = method_label(method)), size = 6.5 * pt, hjust = 1, colour = "grey10") +
  # chip of the colour of the method, before its name
  ggforce::geom_shape(data = bind_rows(lapply(seq_len(nrow(lab_rows)), function(k)
    data.frame(id = k, px = -1.55 + c(0, 0.13, 0.13, 0), py = lab_rows$y[k] + c(-0.22, -0.22, 0.22, 0.22), m = lab_rows$method[k]))),
    aes(x = px, y = py, group = id), fill = rep(unname(method_colors[lab_rows$method]), each = 4), radius = unit(0.6, "pt")) +
  scale_fill_gradientn(colours = heat_v2_cols, limits = c(-lim_v2, lim_v2), oob = scales::squish, breaks = c(-0.04, -0.02, 0, 0.02, 0.04),
                       labels = c("-.04", "-.02", "0", "+.02", "+.04"), name = "Score difference (saturates at +/-.04)",
                       guide = guide_colourbar(barwidth = unit(3.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5,
                                               frame.colour = NA, ticks = FALSE)) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35), guide = "none") +
  scale_colour_identity() +
  scale_x_continuous(limits = c(-2.1, max(col_x) + 0.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(row_y(int_methods)) - 0.55, 0.6), expand = c(0, 0)) +
  labs(title = "Peak: best integrated - best onlyDNAm") +
  theme_void(base_size = 7) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5, margin = margin(b = 4)), legend.position = "bottom",
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0))
figure5_heatmap_v2 = wrap_plots(wrap_plots(p_A + labs(tag = "A") + tag5, p_B2 + labs(tag = "B") + tag5, nrow = 1, widths = c(0.62, 1)),
                                p_C + labs(tag = "C") + tag5, ncol = 1, heights = c(0.8, 1.25))

# ---- supp_generalist_all -------------------------------------------------------------------------------------------
mod_ei = c("preprocessing_mixRNA", "feature_selection_mixRNA", "preprocessing_RNA", "feature_selection_RNA",
           "preprocessing_mixMET", "feature_selection_mixMET", "preprocessing_MET", "feature_selection_MET", "deconvolution")
mod_li = c(setdiff(mod_ei, "deconvolution"), "deconvolution_rna", "deconvolution_met")
gen_all = bind_rows(
  ei_ds %>% mutate(method = as.character(early_integration), id = do.call(paste, c(lapply(across(all_of(mod_ei)), as.character), sep = "|"))),
  li_ds %>% mutate(method = as.character(late_integration), id = do.call(paste, c(lapply(across(all_of(mod_li)), as.character), sep = "|")))) %>%
  select(method, id, dataset, score) %>% mutate(dataset = as.character(dataset)) %>% distinct(method, id, dataset, .keep_all = TRUE)
ds_bench = setdiff(dataset_levels, "VIVO")
u_best = gen_all %>% filter(method %in% c("OnlyMet", "OnlyRna"), dataset %in% ds_bench) %>% group_by(dataset) %>% summarise(u = max(score), .groups = "drop")
gen_sel = gen_all %>% filter(dataset %in% ds_bench) %>% inner_join(u_best, by = "dataset") %>%
  group_by(method, id) %>% filter(n_distinct(dataset) == length(ds_bench)) %>%
  summarise(mean_gain = mean((score - u) / u), .groups = "drop") %>%
  group_by(method) %>% slice_max(mean_gain, n = 1, with_ties = FALSE) %>% ungroup()
# check against the generalist pipelines of 11_exploration_nonlinear_deconv.Rmd (same criterion)
if (file.exists("11_results/top1_recipes.csv")) { # not in script_figure/, which reads 00_results only
  rec11 = read.csv("11_results/top1_recipes.csv")
  chk = gen_sel %>% inner_join(rec11 %>% select(method, mean_gain_11 = mean_gain), by = "method")
  cat(sprintf("check: mean gain of the generalist pipelines vs 11_results/top1_recipes.csv, max abs difference = %.2e\n",
              max(abs(chk$mean_gain - chk$mean_gain_11))))
}
gen_scores = gen_all %>% semi_join(gen_sel, by = c("method", "id"))
ref_dnam = gen_scores %>% filter(method == "OnlyMet") %>% select(dataset, ref = score)
gen_rows = c(int_methods, "OnlyRna")
dG = gen_scores %>% filter(method %in% gen_rows) %>% inner_join(ref_dnam, by = "dataset") %>%
  mutate(diff = score - ref, dataset = factor(dataset, levels = dataset_levels)) %>%
  group_by(dataset) %>% mutate(best_of_dataset = method %in% int_methods & diff == max(diff[method %in% int_methods])) %>% ungroup()

row_yG = function(m) { i = match(as.character(m), gen_rows)
  -(i + ifelse(i > length(early_methods), 0.45, 0) + ifelse(i > length(int_methods), 0.45, 0)) } # gaps: early / late / uni-omic
cellsG = bind_rows(
  dG %>% transmute(method = as.character(method), col = as.character(dataset), v = diff, win = best_of_dataset & dataset != "VIVO"),
  dG %>% filter(dataset != "VIVO") %>% group_by(method) %>% summarise(v = median(diff), .groups = "drop") %>%
    transmute(method = as.character(method), col = "median", v, win = FALSE)) %>%
  mutate(x = col_x[col], y = row_yG(method), faded = col == "VIVO", bold = win | col == "median",
         ink = ifelse(faded, "grey35", ifelse(abs(v) >= 0.7 * lim_v2, "white", "grey15")))
polyG = rect_poly(cellsG); polyG_win = rect_poly(filter(cellsG, win))
lab_rowsG = data.frame(method = gen_rows) %>% mutate(y = row_yG(method))
groupsG = data.frame(y = c(row_yG(early_methods[1]), row_yG(late_methods[1]), row_yG("OnlyRna")) + 0.74,
                     label = c("Early integration", "Late integration", "Uni-omic"),
                     colour = c(unname(method_colors[c(early_methods[4], late_methods[3])]), "grey35"))
chip_col = ifelse(lab_rowsG$method %in% names(method_colors), method_colors[lab_rowsG$method], "grey70")
p_G = ggplot() +
  ggforce::geom_shape(data = polyG, aes(x = px, y = py, group = id, fill = v, alpha = faded), radius = unit(1.2, "pt")) +
  ggforce::geom_shape(data = polyG_win, aes(x = px, y = py, group = id), fill = NA, colour = "grey10", linewidth = 0.5, radius = unit(1.2, "pt")) +
  geom_text(data = cellsG, aes(x = x, y = y, label = fmt_v2(v), colour = ink, fontface = ifelse(bold, "bold", "plain")), size = 6 * pt) +
  geom_text(data = heads, aes(x = x, y = 0.12, label = label, fontface = face, colour = colour), size = 6.3 * pt, vjust = 0) +
  geom_text(data = groupsG, aes(x = -2.05, y = y, label = label, colour = colour), size = 6.5 * pt, fontface = "bold", hjust = 0) +
  geom_text(data = lab_rowsG, aes(x = x_name, y = y, label = method_label(method)), size = 6.5 * pt, hjust = 1, colour = "grey10") +
  ggforce::geom_shape(data = bind_rows(lapply(seq_len(nrow(lab_rowsG)), function(k)
    data.frame(id = k, px = -1.55 + c(0, 0.13, 0.13, 0), py = lab_rowsG$y[k] + c(-0.22, -0.22, 0.22, 0.22)))),
    aes(x = px, y = py, group = id), fill = rep(unname(chip_col), each = 4), radius = unit(0.6, "pt")) +
  scale_fill_gradientn(colours = heat_v2_cols, limits = c(-lim_v2, lim_v2), oob = scales::squish, breaks = c(-0.04, -0.02, 0, 0.02, 0.04),
                       labels = c("-.04", "-.02", "0", "+.02", "+.04"), name = "Score difference (saturates at +/-.04)",
                       guide = guide_colourbar(barwidth = unit(3.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5,
                                               frame.colour = NA, ticks = FALSE)) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35), guide = "none") +
  scale_colour_identity() +
  scale_x_continuous(limits = c(-2.1, max(col_x) + 0.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(row_yG(gen_rows)) - 0.55, 0.6), expand = c(0, 0)) +
  labs(title = "Generalist pipeline - generalist onlyDNAm pipeline") +
  theme_void(base_size = 7) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5, margin = margin(b = 4)), legend.position = "bottom",
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0))
ggsave(file.path(fig_dir, "SuppFig11.pdf"), plot = p_G, width = 13, height = 9, units = "cm", device = cairo_pdf)
