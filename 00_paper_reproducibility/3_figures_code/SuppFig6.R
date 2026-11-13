# Supplementary Figure 6: combinations of RNA and DNAm options within each integration strategy.


data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_02_figure4.R")

# ---- heatmaps_methods_data -----------------------------------------------------------------------------------------
module_method_long = function(df, method_col, modules, methods, reference) {
  df = droplevels(df)
  regime_med = median(df$score_aggreg, na.rm = TRUE)
  out = do.call(rbind, lapply(methods, function(m) {
    dm = df[df[[method_col]] == m, ]
    row_med = median(dm$score_aggreg, na.rm = TRUE)
    ref = if (reference == "row") row_med else regime_med
    do.call(rbind, lapply(names(modules), function(mod) {
      meds = tapply(dm$score_aggreg, droplevels(dm[[modules[[mod]]]]), median, na.rm = TRUE)
      data.frame(method = m, module = mod, level = names(meds), median = as.numeric(meds),
                 value = as.numeric(meds) - ref, row_median = row_med, regime_median = regime_med,
                 stringsAsFactors = FALSE)
    }))
  }))
  out$module = factor(out$module, levels = names(modules))
  out$level = factor(out$level, levels = level_levels)
  out$method = factor(out$method, levels = methods)
  out
}


long_early = list(regime = module_method_long(ei, "early_integration", modules_early, early_methods, "regime"),
                  row = module_method_long(ei, "early_integration", modules_early, early_methods, "row"))
long_late = list(regime = module_method_long(li, "late_integration", modules_late, late_methods, "regime"),
                 row = module_method_long(li, "late_integration", modules_late, late_methods, "row"))

offsets = rbind(
  cbind(regime = "early", unique(long_early$row[, c("method", "row_median", "regime_median")])),
  cbind(regime = "late", unique(long_late$row[, c("method", "row_median", "regime_median")]))
)
offsets$offset_row_minus_regime = round(offsets$row_median - offsets$regime_median, 3)

# ---- heatmaps_methods_fun ------------------------------------------------------------------------------------------
ref_subtitle = c(regime = "difference to the median of all pipelines of the regime",
                 row = "difference to the median of the pipelines of the same method")


# text colour readable on a fill colour (white on the dark fills)
text_on = function(fill) { rgb = grDevices::col2rgb(fill) / 255; ifelse(0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ] < 0.5, "white", "black") }
# cell_h, cell_w (cm): fixed size of the cells (figure 4: the same in B and C); NULL: the cells fill the plot
heatmap_methods = function(long, title, subtitle, methods, modules, lim, sig = NULL, cell_h = NULL, cell_w = NULL, sig_star = FALSE) {
  strip_cols = unname(module_pal[module_family[names(modules)]])
  meth_cols = unname(method_colors[methods])
  # columns: one strip per module (colour of the module); rows: one strip per method (colour of the method), on the left
  facet_layer = ggh4x::facet_grid2(method ~ module, scales = "free", space = "free", switch = "y",
                                   labeller = labeller(module = function(x) sub(" ", "\n", x), method = function(x) method_label(x)),
                                   strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_cols, colour = NA),
                                                               background_y = ggh4x::elem_list_rect(fill = meth_cols, colour = NA),
                                                               text_y = ggh4x::elem_list_text(colour = text_on(meth_cols))))
  ggplot(long, aes(x = level, y = method, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_layer +
    (if (!is.null(sig) && !sig_star) geom_point(data = filter(sig, sig), aes(x = level, y = method), size = 0.7, colour = "black", inherit.aes = FALSE)) +
    # sig_star = TRUE: a star rather than a dot (no confusion with the circles of the interactions)
    (if (!is.null(sig) && sig_star) geom_text(data = filter(sig, sig), aes(x = level, y = method), label = "*", size = 10 / .pt, vjust = 0.75, fontface = "bold",
                                              colour = "grey10", inherit.aes = FALSE)) +
    (if (!is.null(cell_h)) ggh4x::force_panelsizes(rows = unit(rep(cell_h, length(methods)), "cm"),
                                                   cols = unit(cell_w * sapply(names(modules), function(mo) n_distinct(long$level[long$module == mo])), "cm"))) +
    fill_scale(lim) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey60", fill = NA, linewidth = 0.3),
      panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(1, "pt"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(size = 6.5, face = "bold", colour = "black", lineheight = 0.9, margin = margin(2, 1, 2, 1)),
      strip.text.y.left = element_text(size = 6.5, face = "bold", angle = 0, hjust = 1, margin = margin(1, 3, 1, 3)),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      plot.subtitle = element_text(hjust = 0.5, size = 7, colour = "grey30")
    )
}

# ---- panels_CH_fun -------------------------------------------------------------------------------------------------
module_violin = function(data, col, title, show_y = TRUE, base_size = 8, adjust = 1.3, score = "score_aggreg", letters_lv = NULL, level_order = NULL) {
  data[[col]] = droplevels(data[[col]])
  med_by_level = tapply(data[[score]], data[[col]], median, na.rm = TRUE)
  lv = if (!is.null(level_order)) level_order else names(sort(med_by_level, decreasing = TRUE))
  data[[col]] = factor(as.character(data[[col]]), levels = lv)
  sm = data.frame(level = factor(lv, levels = lv), med = as.numeric(med_by_level[lv]))

  ggplot(data, aes(x = .data[[col]], y = .data[[score]], fill = .data[[col]])) +
    geom_violin(trim = TRUE, scale = "width", adjust = adjust, width = 0.9, colour = "grey30", linewidth = 0.25) +
    geom_errorbar(data = sm, aes(x = level, ymin = med, ymax = med), width = 0.6, linewidth = 0.6,
                  colour = "black", inherit.aes = FALSE) +
    geom_text(data = sm, aes(x = level, y = 1.01, label = sprintf("%.2f", med)), size = base_size / 4.5, angle = 90,
              hjust = 0, inherit.aes = FALSE) + # median written above each violin (vertical: no overlap)
    (if (!is.null(letters_lv)) geom_text(data = sm %>% mutate(l = letters_lv[as.character(level)]), aes(x = level, y = -0.07, label = l),
                                         size = base_size / 3.2, fontface = "bold", colour = "grey15", inherit.aes = FALSE)) +
    scale_fill_manual(values = level_colors[lv]) +
    scale_y_continuous(limits = c(-0.1, 1.15), breaks = seq(0, 1, 0.25)) + # letters of the tests under the violins
    guides(fill = "none") +
    labs(title = title, x = NULL, y = "Aggregated score") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.4),
      axis.line = element_line(colour = "grey40", linewidth = 0.4),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size - 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size)
    ) +
    (if (show_y) NULL else theme(axis.title.y = element_blank(), axis.text.y = element_blank()))
}

interaction_tile = function(df, rna_col, met_col, rna_lab, met_lab, lim = lim_modules, sig = NULL) {
  df = droplevels(df)
  baseline = median(df$score_aggreg, na.rm = TRUE)
  m = tapply(df$score_aggreg, list(df[[rna_col]], df[[met_col]]), median, na.rm = TRUE) - baseline
  long = data.frame(rna = rep(rownames(m), times = ncol(m)), met = rep(colnames(m), each = nrow(m)), value = as.vector(m),
                    stringsAsFactors = FALSE)
  long$rna = factor(long$rna, levels = rev(rownames(m)))
  long$met = factor(long$met, levels = colnames(m))
  if (!is.null(sig)) long = long %>% left_join(sig %>% select(rna, met, sig), by = c("rna", "met")) %>%
    mutate(rna = factor(rna, levels = rev(rownames(m))), met = factor(met, levels = colnames(m)))
  ggplot(long, aes(x = met, y = rna, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    (if (!is.null(sig)) geom_point(data = function(x) filter(x, sig %in% TRUE), size = 1, colour = "black")) +
    fill_scale(lim) +
    labs(x = met_lab, y = rna_lab) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title = element_text(size = 7, face = "bold")
    )
}

# heading of a row of panels: coloured band (module family) with the name of the module
module_heading = function(txt, family) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = module_pal[[family]]) +
    annotate("text", x = 0.5, y = 0.5, label = txt, fontface = "bold", size = 3.4, colour = "black") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    theme_void()
}
heading_height = 0.13 # height of the coloured bands relative to the panels below
n_lev = function(df, col) nlevels(droplevels(df[[col]]))

cell_modules = function(df, rna_col, met_col, rna_lab, met_lab, heading, family, tag, regime = NULL) {
  nR = n_lev(df, rna_col); nM = n_lev(df, met_col)
  lt = function(cl) if (!is.null(regime)) mstats[[regime]][[cl]]$letters
  body = wrap_plots(
    module_violin(df, rna_col, "RNA", show_y = TRUE, letters_lv = lt(rna_col)) + labs(tag = tag) + tag_theme,
    module_violin(df, met_col, "DNAm", show_y = FALSE, letters_lv = lt(met_col)),
    interaction_tile(df, rna_col, met_col, rna_lab, met_lab, sig = if (!is.null(regime)) istats[[regime]][[paste(rna_col, met_col, sep = "|")]]),
    nrow = 1, widths = c(nR + 1.6, nM + 0.3, 0.8 * nM + 2.6)
  )
  wrap_plots(module_heading(heading, family), body, ncol = 1, heights = c(heading_height, 1))
}
cell_deconv_early = function(df, heading, family, tag) {
  n = n_lev(df, "deconvolution")
  body = wrap_plots(
    module_violin(df, "deconvolution", "RNA + DNAm (joint)", show_y = TRUE, letters_lv = mstats$early$deconvolution$letters) + labs(tag = tag) + tag_theme,
    plot_spacer(), plot_spacer(), nrow = 1, widths = c(n + 1.6, n + 0.3, 0.8 * n + 2.6)
  )
  wrap_plots(module_heading(heading, family), body, ncol = 1, heights = c(heading_height, 1))
}

# ---- stats ---------------------------------------------------------------------------------------------------------
ds_tab = list(early = read_res(file.path(data_dir, "ei_ds.rds")) %>% filter(method %in% early_methods, dataset != "VIVO") %>% droplevels(),
              late = read_res(file.path(data_dir, "li_ds.rds")) %>% filter(method %in% late_methods, dataset != "VIVO") %>% droplevels())
key_all = list(early = all_functions_ei, late = all_functions_li)
mix_of = c(preprocessing_RNA = "preprocessing_mixRNA", feature_selection_RNA = "feature_selection_mixRNA",
           preprocessing_MET = "preprocessing_mixMET", feature_selection_MET = "feature_selection_mixMET")
paste_key = function(d, cols) do.call(paste, c(lapply(d[cols], as.character), sep = "|"))
# exact two-sided Wilcoxon signed-rank p-value, by enumeration of the 2^n sign permutations: zero values are dropped and
# ties get mid-ranks (wilcox.test(exact = TRUE) falls back to the normal approximation as soon as there are ties or zeros)
exact_signrank_p = function(x) {
  x = x[x != 0]; n = length(x)
  if (n == 0) return(1)
  r = rank(abs(x)); t_obs = sum(r[x > 0]); e = sum(r) / 2
  signs = as.matrix(expand.grid(rep(list(c(0, 1)), n)))
  t_all = drop(signs %*% r)
  mean(abs(t_all - e) >= abs(t_obs - e) - 1e-9)
}
w_test = function(x) if (length(x) >= 5 && any(x != 0)) exact_signrank_p(x) else 1

# letters of a compact letter display, relabelled so that the level with the highest median gets "a"
relabel_letters = function(lett, ord) {
  lett = lett[ord]; old = unique(unlist(strsplit(paste(lett, collapse = ""), "")))
  map = setNames(letters[seq_along(old)], old)
  sapply(lett, function(x) paste(sort(map[strsplit(x, "")[[1]]]), collapse = ""))
}
module_stats = function(regime, col) {
  d = ds_tab[[regime]]
  keys = setdiff(key_all[[regime]], c(col, mix_of[col]))
  d = d %>% mutate(key = paste_key(., keys), level = as.character(.data[[col]]))
  lev_med = d %>% group_by(dataset, level) %>% summarise(m = median(score), .groups = "drop") %>% tidyr::pivot_wider(names_from = level, values_from = m)
  lv = setdiff(names(lev_med), "dataset"); mat = as.matrix(lev_med[, lv])
  fr = friedman.test(mat)
  pw = bind_rows(lapply(combn(lv, 2, simplify = FALSE), function(pr) {
    md = inner_join(d %>% filter(level == pr[1]) %>% select(dataset, key, sa = score),
                    d %>% filter(level == pr[2]) %>% select(dataset, key, sb = score), by = c("dataset", "key")) %>%
      group_by(dataset) %>% summarise(diff = median(sa - sb), n = n(), .groups = "drop")
    data.frame(level_a = pr[1], level_b = pr[2], matched_pairs_per_dataset = median(md$n), median_diff_a_minus_b = median(md$diff),
               n_datasets_a_better = sum(md$diff > 0), p = w_test(md$diff))
  })) %>% mutate(p_BH = p.adjust(p, method = "BH"))
  # order of the levels from the matched comparisons: sum of the median matched differences against the other levels
  lv_score = sapply(lv, function(l) sum(pw$median_diff_a_minus_b[pw$level_a == l]) - sum(pw$median_diff_a_minus_b[pw$level_b == l]))
  ord = names(sort(lv_score, decreasing = TRUE))
  lett = multcompView::multcompLetters(setNames(pw$p_BH < 0.05, paste(pw$level_a, pw$level_b, sep = "-")))$Letters
  list(friedman = data.frame(regime = regime, module = col, n_levels = length(lv), friedman_chi2 = unname(fr$statistic), p = fr$p.value,
                             kendall_W = unname(fr$statistic) / (nrow(mat) * (length(lv) - 1))),
       pairwise = pw %>% mutate(regime = regime, module = col, .before = 1),
       letters = relabel_letters(lett, ord), order = ord)
}
interaction_stats = function(regime, rna_col, met_col) {
  per = ds_tab[[regime]] %>% mutate(rna = as.character(.data[[rna_col]]), met = as.character(.data[[met_col]])) %>%
    group_by(dataset) %>% mutate(g = median(score)) %>%
    group_by(dataset, rna) %>% mutate(r = median(score)) %>% group_by(dataset, met) %>% mutate(cc = median(score)) %>%
    group_by(dataset, rna, met) %>% summarise(contrast = median(score) - first(r) - first(cc) + first(g), .groups = "drop")
  per %>% group_by(rna, met) %>% summarise(interaction = median(contrast), p = w_test(contrast), .groups = "drop") %>%
    mutate(p_BH = p.adjust(p, method = "BH"), sig = p_BH < 0.05)
}
method_cell_stats = function(regime, methods, modules) {
  d = ds_tab[[regime]]
  out = bind_rows(lapply(methods, function(m) bind_rows(lapply(names(modules), function(mod) {
    per = d %>% filter(method == m) %>% group_by(dataset) %>% mutate(ref = median(score)) %>%
      group_by(dataset, level = as.character(.data[[modules[[mod]]]])) %>% summarise(dev = median(score) - first(ref), .groups = "drop")
    wide = per %>% tidyr::pivot_wider(names_from = level, values_from = dev)
    fr = friedman.test(as.matrix(wide[, -1]))
    per %>% group_by(level) %>% summarise(median_dev = median(dev), p = w_test(dev), .groups = "drop") %>%
      mutate(method = m, module = mod, kendall_W = unname(fr$statistic) / (nrow(wide) * (ncol(wide) - 2)), friedman_p = fr$p.value)
  }))))
  out %>% mutate(p_BH = p.adjust(p, method = "BH"), sig = p_BH < 0.05,
                 module = factor(module, levels = names(modules)), level = factor(level, levels = level_levels),
                 method = factor(method, levels = methods))
}

mod_cols = list(early = unname(unlist(modules_early)), late = unname(unlist(modules_late)))
mstats = lapply(c(early = "early", late = "late"), function(rg) setNames(lapply(mod_cols[[rg]], function(cl) module_stats(rg, cl)), mod_cols[[rg]]))
int_pairs = list(early = list(c("preprocessing_RNA", "preprocessing_MET"), c("feature_selection_RNA", "feature_selection_MET")),
                 late = list(c("preprocessing_RNA", "preprocessing_MET"), c("feature_selection_RNA", "feature_selection_MET"), c("deconvolution_rna", "deconvolution_met")))
istats = lapply(c(early = "early", late = "late"), function(rg) setNames(lapply(int_pairs[[rg]], function(pr) interaction_stats(rg, pr[1], pr[2])),
                                                                        sapply(int_pairs[[rg]], paste, collapse = "|")))
cstats = list(early = method_cell_stats("early", early_methods, modules_early), late = method_cell_stats("late", late_methods, modules_late))

friedman_tab = bind_rows(lapply(mstats, function(x) bind_rows(lapply(x, `[[`, "friedman"))))
pairwise_tab = bind_rows(lapply(mstats, function(x) bind_rows(lapply(x, `[[`, "pairwise"))))

# ---- variance ------------------------------------------------------------------------------------------------------
var_decomp = function(regime, B = 100) {
  d = ds_tab[[regime]] %>% mutate(dataset = factor(dataset))
  terms = c("method", mod_cols[[regime]])
  f = as.formula(paste("score ~ dataset + (", paste(terms, collapse = " + "), ")^2"))
  eta = function(dd) { a = anova(lm(f, data = dd)); ss = setNames(a[["Sum Sq"]], rownames(a)); ss[setdiff(names(ss), c("dataset", "Residuals"))] / sum(ss) }
  est = eta(d)
  set.seed(1)
  boot = sapply(seq_len(B), function(b) {
    pick = sample(levels(d$dataset), replace = TRUE)
    dd = bind_rows(lapply(seq_along(pick), function(i) d %>% filter(dataset == pick[i]) %>% mutate(dataset = paste0(pick[i], "_", i)))) %>% mutate(dataset = factor(dataset))
    e = eta(dd); e[names(est)]
  })
  data.frame(regime = regime, term = names(est), eta2 = unname(est), low = apply(boot, 1, quantile, 0.025, na.rm = TRUE), high = apply(boot, 1, quantile, 0.975, na.rm = TRUE))
}
short = c(method = "integration strategy", preprocessing_RNA = "PP RNA", preprocessing_MET = "PP DNAm", feature_selection_RNA = "FS RNA",
          feature_selection_MET = "FS DNAm", deconvolution = "Deconv", deconvolution_rna = "Deconv RNA", deconvolution_met = "Deconv DNAm")

# ---- options_data --------------------------------------------------------------------------------------------------
# per-dataset median of the matched differences (level - reference level) for one module
matched_vs_ref = function(regime, col) {
  d = ds_tab[[regime]]
  keys = setdiff(key_all[[regime]], c(col, mix_of[col]))
  d = d %>% mutate(key = paste_key(., keys), level = as.character(.data[[col]]))
  ord = mstats[[regime]][[col]]$order; ref_level = ord[1] # order from the matched comparisons; best option = reference (not "ref": the tables have a column ref)
  r = d %>% filter(level == ref_level) %>% select(dataset, key, s_ref = score)
  bind_rows(lapply(ord, function(l) d %>% filter(level == l) %>% select(dataset, key, s = score) %>% inner_join(r, by = c("dataset", "key")) %>%
                     group_by(dataset) %>% summarise(diff = median(s - s_ref), .groups = "drop") %>% mutate(level = l))) %>%
    mutate(level = factor(level, levels = ord), reference = ref_level)
}
mdiff = lapply(c(early = "early", late = "late"), function(rg) setNames(lapply(mod_cols[[rg]], function(cl) matched_vs_ref(rg, cl)), mod_cols[[rg]]))

paired_plot = function(regime, col, title, show_y = TRUE, base_size = 8, show_title = TRUE) {
  d = mdiff[[regime]][[col]]; lt = mstats[[regime]][[col]]$letters
  med = d %>% group_by(level) %>% summarise(m = median(diff), .groups = "drop") %>% mutate(l = lt[as.character(level)])
  rng = range(c(d$diff, 0)); pad = 0.12 * diff(rng) + 0.002
  ggplot(d, aes(level, diff)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_point(aes(fill = level), shape = 21, size = 1.3, stroke = 0.2, position = position_jitter(width = 0.12, height = 0, seed = 1)) +
    geom_errorbar(data = med, aes(x = level, ymin = m, ymax = m), width = 0.55, linewidth = 0.6, inherit.aes = FALSE) +
    geom_text(data = med, aes(x = level, y = rng[1] - pad, label = l), size = base_size / 3.2, fontface = "bold", inherit.aes = FALSE) +
    scale_fill_manual(values = level_colors[levels(d$level)], guide = "none") +
    scale_y_continuous(limits = c(rng[1] - 1.6 * pad, rng[2] + 0.4 * pad)) +
    scale_x_discrete(expand = expansion(add = 0.6)) +
    labs(title = if (show_title) title, x = NULL, y = "Change vs\nbest option") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
          axis.line = element_line(colour = "grey40", linewidth = 0.4), axis.title.y = element_text(size = base_size - 1.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size - 1), plot.title = element_text(hjust = 0.5, face = "bold", size = base_size)) +
    (if (show_y) NULL else theme(axis.title.y = element_blank()))
}

stars = function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
# pairwise tests: for two options of a module (columns of the supplementary figure of the pairwise tests),
# one point per dataset = median of the score differences between matched pipelines (identical for all the other
# modules): exactly the 8 values of the pairwise test; stars = BH-adjusted p of this pair (same test as the letters).
pair_dataset_diffs = function(regime, col, a, b) {
  keys = setdiff(key_all[[regime]], c(col, mix_of[col]))
  d = ds_tab[[regime]] %>% mutate(key = paste_key(., keys), level = as.character(.data[[col]]))
  inner_join(d %>% filter(level == a) %>% select(dataset, key, sa = score), d %>% filter(level == b) %>% select(dataset, key, sb = score), by = c("dataset", "key")) %>%
    group_by(dataset) %>% summarise(v = median(sa - sb), .groups = "drop")
}
violin_ordered_top = function(df, regime, col, title, show_y) module_violin(df, col, title, show_y = show_y, level_order = mstats[[regime]][[col]]$order) +
  scale_x_discrete(expand = expansion(add = 0.6)) + labs(y = "Aggregated\nscore") + theme(axis.text.x = element_blank())

# one row of panels (RNA, DNAm, interaction heatmap with the dots of the interaction tests), violins above the matched differences
cell_option = function(df, regime, rna_col, met_col, rna_lab, met_lab, heading, family, tag) {
  nR = n_lev(df, rna_col); nM = n_lev(df, met_col)
  heat = interaction_tile(df, rna_col, met_col, rna_lab, met_lab, sig = istats[[regime]][[paste(rna_col, met_col, sep = "|")]])
  w = c(nR + 1.6, nM + 0.3, 0.8 * nM + 2.6)
  body = wrap_plots(violin_ordered_top(df, regime, rna_col, "RNA", TRUE) + labs(tag = tag) + tag_theme, violin_ordered_top(df, regime, met_col, "DNAm", FALSE), heat,
                    paired_plot(regime, rna_col, "", TRUE, show_title = FALSE), paired_plot(regime, met_col, "", FALSE, show_title = FALSE), plot_spacer(),
                    design = "ABC\nDEC", widths = w, heights = c(1, 0.75))
  wrap_plots(module_heading(heading, family), body, ncol = 1, heights = c(heading_height, 1))
}
cell_option_deconv_early = function(heading, family, tag) {
  n = n_lev(ei, "deconvolution")
  w = c(n + 1.6, n + 0.3, 0.8 * n + 2.6)
  body = wrap_plots(violin_ordered_top(ei, "early", "deconvolution", "RNA + DNAm (joint)", TRUE) + labs(tag = tag) + tag_theme, plot_spacer(), plot_spacer(),
                    paired_plot("early", "deconvolution", "", show_title = FALSE), plot_spacer(), plot_spacer(), nrow = 2, widths = w, heights = c(1, 0.75))
  wrap_plots(module_heading(heading, family), body, ncol = 1, heights = c(heading_height, 1))
}
compose_violins = function() {
  top = wrap_plots(
    heatmap_methods(long_early$row, "Early integration", ref_subtitle[["row"]], early_methods, modules_early, lim_top[["row"]], sig = cstats$early) + labs(tag = "A") + tag_theme,
    heatmap_methods(long_late$row, "Late integration", ref_subtitle[["row"]], late_methods, modules_late, lim_top[["row"]], sig = cstats$late) + labs(tag = "B") + tag_theme,
    nrow = 1, widths = c(23, 30))
  row = function(df, rg, rc, mc, rl, ml, h, fam, tg) cell_option(df, rg, rc, mc, rl, ml, h, fam, tg)
  row_pp = wrap_plots(row(ei, "early", "preprocessing_RNA", "preprocessing_MET", "Preprocessing RNA", "Preprocessing DNAm", "Early integration - Pre-processing", "preprocessing", "C"),
                      row(li, "late", "preprocessing_RNA", "preprocessing_MET", "Preprocessing RNA", "Preprocessing DNAm", "Late integration - Pre-processing", "preprocessing", "D"), nrow = 1)
  row_fs = wrap_plots(row(ei, "early", "feature_selection_RNA", "feature_selection_MET", "FS RNA", "FS DNAm", "Early integration - Feature selection", "feature_selection", "E"),
                      row(li, "late", "feature_selection_RNA", "feature_selection_MET", "FS RNA", "FS DNAm", "Late integration - Feature selection", "feature_selection", "F"), nrow = 1)
  row_dc = wrap_plots(cell_option_deconv_early("Early integration - Deconvolution", "deconvolution", "G"),
                      row(li, "late", "deconvolution_rna", "deconvolution_met", "Deconvolution RNA", "Deconvolution DNAm", "Late integration - Deconvolution", "deconvolution", "H"), nrow = 1)
  wrap_plots(top, row_pp, row_fs, row_dc, ncol = 1, heights = c(1.05, 2.2, 2.3, 2.5), guides = "collect") & theme(legend.position = "bottom")
}

# ---- option_E_data -------------------------------------------------------------------------------------------------
interaction_stats_method = function(regime, m, rna_col, met_col) {
  ds_tab[[regime]] %>% filter(method == m) %>%
    mutate(rna = as.character(.data[[rna_col]]), met = as.character(.data[[met_col]])) %>%
    # mean-based interaction contrast (exact two-way decomposition; the median version is not additive with bimodal scores)
    group_by(dataset) %>% mutate(g = mean(score)) %>%
    group_by(dataset, rna) %>% mutate(r = mean(score)) %>% group_by(dataset, met) %>% mutate(cc = mean(score)) %>%
    group_by(dataset, rna, met) %>% summarise(contrast = mean(score) - first(r) - first(cc) + first(g), .groups = "drop") %>%
    group_by(rna, met) %>% summarise(interaction = median(contrast), p = w_test(contrast), .groups = "drop") %>%
    mutate(p_BH = p.adjust(p, method = "BH"), sig = p_BH < 0.05)
}
pair_defs = list(early = list("Pre-processing" = c("preprocessing_RNA", "preprocessing_MET"),
                              "Feature selection" = c("feature_selection_RNA", "feature_selection_MET")),
                 late = list("Pre-processing" = c("preprocessing_RNA", "preprocessing_MET"),
                             "Feature selection" = c("feature_selection_RNA", "feature_selection_MET"),
                             "Deconvolution" = c("deconvolution_rna", "deconvolution_met")))
method_interactions = function(regime) {
  df = if (regime == "early") ei else li
  mcol = if (regime == "early") "early_integration" else "late_integration"
  methods = if (regime == "early") early_methods else late_methods
  bind_rows(lapply(methods, function(m) bind_rows(lapply(names(pair_defs[[regime]]), function(pn) {
    rc = pair_defs[[regime]][[pn]][1]; mc = pair_defs[[regime]][[pn]][2]
    dm = droplevels(df[df[[mcol]] == m, ])
    tab = tapply(dm$score_aggreg, list(as.character(dm[[rc]]), as.character(dm[[mc]])), median, na.rm = TRUE) - median(dm$score_aggreg, na.rm = TRUE)
    long = data.frame(rna = rep(rownames(tab), times = ncol(tab)), met = rep(colnames(tab), each = nrow(tab)), value = as.vector(tab))
    long %>% left_join(interaction_stats_method(regime, m, rc, mc) %>% select(rna, met, interaction, p_BH, sig), by = c("rna", "met")) %>%
      mutate(method = m, pair = pn, rna_order = match(rna, level_levels), met_order = match(met, level_levels))
  }))))
}
mint = lapply(c(early = "early", late = "late"), method_interactions)

# ---- figure4 -------------------------------------------------------------------------------------------------------
# drawn at its final size (18 cm wide, portrait A4 page): the font sizes are the printed sizes, at least 6 pt
tagE = theme(plot.tag = element_text(face = "bold", size = 14), plot.background = element_blank()) # transparent: the larger tags must not hide the axis title of A
lim_int = 0.03 # colour limit of the interaction contrasts (D, E)
cell_cm = 0.24 # size of a cell of the interaction heatmaps (cm), identical in D and E
pair_family = c("Pre-processing" = "preprocessing", "Feature selection" = "feature_selection", "Deconvolution" = "deconvolution")
# fill_col: colour of the cells, "interaction" (interaction contrast) or "value" (score of the combination: median of the
# cell - median of the method); lim: limit of the colour scale. bubbles = TRUE (figure 4): in each cell, an empty circle of
# area proportional to |interaction contrast| (scale shared by D and E, up to int_size_max), and "+" / "-" for the
# significant positive / negative interactions; bubbles = FALSE: a dot for the significant interactions
small_heat2 = function(regime, m, pn, first_col, first_row, last_col = FALSE, top_pair = FALSE, bottom_pair = FALSE,
                       fill_col = "interaction", lim = lim_int, bubbles = FALSE, sign_fill = FALSE) {
  short = c("Pre-processing" = "PP", "Feature selection" = "FS", "Deconvolution" = "Deconv")
  d = mint[[regime]] %>% filter(method == m, pair == pn) %>%
    mutate(rna = reorder(rna, -rna_order), met = reorder(met, met_order), pair = short[[pn]], meth = method_label(m))
  mc = unname(method_colors[[m]])
  ggplot(d, aes(met, rna, fill = .data[[fill_col]])) +
    geom_tile(colour = "white", linewidth = 0.4) +
    (if (!bubbles) geom_point(data = filter(d, sig %in% TRUE), size = 0.4, colour = "black")
     else if (sign_fill) list( # significant interactions only: filled circle if positive, open (white) circle if negative
               geom_point(data = filter(d, sig %in% TRUE, interaction > 0), aes(size = abs(interaction)), shape = 16, colour = "grey10"),
               geom_point(data = filter(d, sig %in% TRUE, interaction <= 0), aes(size = abs(interaction)), shape = 21, fill = "white", colour = "grey10", stroke = 0.35),
               scale_size_area(limits = c(0, int_size_max), max_size = 2.6, oob = scales::squish, guide = "none"))
     else list(geom_point(aes(size = abs(interaction)), shape = 1, colour = "grey10", stroke = 0.3),
               geom_text(data = filter(d, sig %in% TRUE), aes(label = ifelse(interaction > 0, "+", "-")), size = 6 / .pt, fontface = "bold", vjust = 0.4),
               scale_size_area(limits = c(0, int_size_max), max_size = 2.6, oob = scales::squish, guide = "none"))) +
    scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish, na.value = "grey90", name = "Interaction contrast",
                         guide = "none") + # legend drawn once below the figure
    scale_y_discrete(position = "right") + # RNA option names on the right of the last column
    # "RNA option" on top of the column of RNA option names (last column, first pair); "DNAm option" under the first column (last pair)
    # (a narrow panel is centred in its column, whose width is set by the pair with the most DNAm options: shift to the labels)
    (if (last_col && top_pair) annotate("text", x = nlevels(droplevels(d$met)) + 0.65 + (max(table(distinct(mint[[regime]], pair, met)$pair)) - nlevels(droplevels(d$met))) / 2,
                                        y = nlevels(droplevels(d$rna)) + 0.6, label = "RNA option",
                                        hjust = 0, vjust = 0, size = 6 / .pt, fontface = "bold.italic")) +
    # limits fixed to the cells: the annotation outside the panel does not widen the axes
    coord_cartesian(xlim = c(0.5, nlevels(droplevels(d$met)) + 0.5), ylim = c(0.5, nlevels(droplevels(d$rna)) + 0.5), expand = FALSE, clip = "off") +
    # "DNAm option" at the far left of the row of DNAm option names, vertical like the module strips (first column, last pair)
    (if (first_col && bottom_pair) annotate("text", x = -1.1, y = -2.8, label = "DNAm option", angle = 90, size = 6 / .pt, fontface = "bold.italic")) +
    labs(x = NULL, y = NULL) +
    # rows: pair of modules (strip on the left, first column only); columns: method (strip on top, first row only)
    ggh4x::facet_grid2(rows = vars(pair), cols = vars(meth), switch = "y",
                       strip = if (first_row) ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = unname(module_pal[[pair_family[[pn]]]]), colour = NA),
                                                                   background_x = ggh4x::elem_list_rect(fill = mc, colour = NA),
                                                                   text_x = ggh4x::elem_list_text(colour = text_on(mc)))
                               else ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = unname(module_pal[[pair_family[[pn]]]]), colour = NA))) +
    ggh4x::force_panelsizes(rows = unit(nlevels(droplevels(d$rna)) * cell_cm, "cm"), cols = unit(nlevels(droplevels(d$met)) * cell_cm, "cm")) +
    theme_minimal(base_size = 6.5) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                    axis.text.y = element_blank(), axis.text.y.right = if (last_col) element_text(size = 6, hjust = 0) else element_blank(),
          strip.placement = "outside",
          strip.text.y.left = if (first_col) element_text(angle = 90, face = "bold", size = 6) else element_blank(),
          strip.background.y = if (first_col) NULL else element_blank(),
          strip.text.x = if (first_row) element_text(face = "bold", size = 6, margin = margin(1.5, 1, 1.5, 1)) else element_blank(),
          strip.background.x = if (first_row) NULL else element_blank(),
          plot.margin = margin(2, 1, 4, 1))
}
# header of a column of D / E: the method name on a band of its colour, as wide as the column (long names on two lines)
method_header = function(m) {
  mc = unname(method_colors[[m]])
  lab = method_label(m)
  ggplot() + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = mc) +
    annotate("text", x = 0.5, y = 0.5, label = lab, size = 6 / .pt, fontface = "bold", lineheight = 0.85, colour = text_on(mc)) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    ggh4x::force_panelsizes(rows = unit(0.5, "cm")) + # fixed height; the width is that of the column
    theme_void() + theme(plot.margin = margin(0, 1, 2, 1))
}
interaction_panel2 = function(regime, title, ...) { # ...: fill_col and lim of small_heat2
  methods = if (regime == "early") early_methods else late_methods
  pairs = names(pair_defs[[regime]])
  heads = lapply(methods, method_header)
  plots = unlist(lapply(seq_along(pairs), function(i) lapply(seq_along(methods), function(j)
    small_heat2(regime, methods[j], pairs[i], j == 1, FALSE, j == length(methods), i == 1, i == length(pairs), ...))), recursive = FALSE) # method names in the header row
  # last row: spacers of free height, so that the heatmaps (fixed sizes) sit at the top of their area, the title at the level of the tag
  spacers = lapply(methods, function(m) ggplot() + theme_void() + ggh4x::force_panelsizes(cols = unit(0.1, "cm"))) # width fixed, height free
  wrap_plots(c(heads, plots, spacers), ncol = length(methods), byrow = TRUE) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 7, margin = margin(b = 5, l = 14)))) # away from the tag and the method headers
}
small_heat_theme = theme(plot.title = element_text(size = 7.5),
                         strip.text.x = element_text(size = 6, face = "bold", lineheight = 0.85, margin = margin(1.5, 0, 1.5, 0)))
# D, E: cells coloured by the score of the combination (scale of B and C), empty circle of area proportional to
# |interaction contrast| (largest circle for the contrasts >= 0.03: only concatnoscale pre-processing reaches 0.08, and
# with it as the maximum the other circles were too small), "+" / "-" for the significant interactions
int_size_max = 0.03
int_early = interaction_panel2("early", "Early integration", fill_col = "value", lim = lim_top[["row"]], bubbles = TRUE)
int_late = interaction_panel2("late", "Late integration", fill_col = "value", lim = lim_top[["row"]], bubbles = TRUE)
int_legend = (ggplot(data.frame(x = 1:4, y = 1, s = c(0.005, 0.01, 0.02, 0.03)), aes(x, y, size = s)) + geom_point(shape = 1, alpha = 0) + # legend only
  scale_size_area(limits = c(0, int_size_max), max_size = 2.6, oob = scales::squish, breaks = c(0.005, 0.01, 0.02, 0.03), labels = c("0.005", "0.01", "0.02", ">= 0.03"),
                  name = "|Interaction|", guide = guide_legend(override.aes = list(alpha = 1, colour = "grey10", stroke = 0.3))) +
  theme_void())
# A on the whole width; below, two columns: B over D (early integration), C over E (late integration). The cells of B
# and C have the same size, so that C (3 methods) is shorter than B (5 methods) and E moves up under C.
legend_theme = theme(legend.position = "bottom", legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
                     legend.key.size = unit(0.3, "cm"), legend.margin = margin(0, 4, 0, 4), legend.box.spacing = unit(2, "pt"))
cell_bc = c(h = 0.42, w = 0.24) # cm
# B and C in frames (fixed cell sizes); the colour scale of B-E drawn once, at the bottom, with the circles of D and E
pB4 = heatmap_methods(long_early$row, "Early integration", NULL, early_methods, modules_early, lim_top[["row"]], sig = cstats$early,
                      cell_h = cell_bc[["h"]], cell_w = cell_bc[["w"]]) + small_heat_theme + theme(legend.position = "none")
pC4 = heatmap_methods(long_late$row, "Late integration", NULL, late_methods, modules_late, lim_top[["row"]], sig = cstats$late,
                      cell_h = cell_bc[["h"]], cell_w = cell_bc[["w"]]) + small_heat_theme + theme(legend.position = "none")
fill_legend = (ggplot(data.frame(x = 1, y = 1, v = 0), aes(x, y, fill = v)) + geom_tile(alpha = 0) + # legend only
  scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                       limits = c(-lim_top[["row"]], lim_top[["row"]]), name = "Score difference (B-E)",
                       guide = guide_colourbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5)) +
  theme_void())
col_early = wrap_plots(wrap_elements(full = pB4) + labs(tag = "B") + tagE, wrap_elements(full = int_early) + labs(tag = "D") + tagE,
                       ncol = 1, heights = c(5.4, 8.2))
col_late = wrap_plots(wrap_elements(full = pC4) + labs(tag = "C") + tagE, wrap_elements(full = int_late) + labs(tag = "E") + tagE,
                      ncol = 1, heights = c(4.6, 9.0)) # C (3 methods) shorter than B (5 methods): E moves up

# ---- figure4_v2 ----------------------------------------------------------------------------------------------------
fs_lab = 6; fs_name = 6.5; fs_title = 7.5 # font sizes (pt)
pt = 1 / ggplot2::.pt # geom_text size for 1 pt
# ---- B, C: stars ----
bc_theme = theme(plot.title = element_text(size = fs_title, face = "bold"), axis.text.x = element_text(size = fs_lab),
                 strip.text.x = element_text(size = fs_lab, face = "bold", lineheight = 0.85, margin = margin(1.5, 0, 1.5, 0)),
                 strip.text.y.left = element_text(size = fs_name, face = "bold", angle = 0, hjust = 1, margin = margin(1, 3, 1, 3)),
                 legend.position = "none", plot.background = element_blank()) # transparent: must not hide the axis title of A
pB5 = heatmap_methods(long_early$row, "Early integration", NULL, early_methods, modules_early, lim_top[["row"]], sig = cstats$early,
                      cell_h = cell_bc[["h"]], cell_w = cell_bc[["w"]], sig_star = TRUE) + bc_theme
pC5 = heatmap_methods(long_late$row, "Late integration", NULL, late_methods, modules_late, lim_top[["row"]], sig = cstats$late,
                      cell_h = cell_bc[["h"]], cell_w = cell_bc[["w"]], sig_star = TRUE) + bc_theme
# ---- D, E: RNA options in rows (left, first method only), DNAm options in columns (bottom of each heatmap) ----
cell_de = 0.26 # cm, identical in D and E
small_heat3 = function(regime, m, pn, first_col) {
  short = c("Pre-processing" = "PP", "Feature selection" = "FS", "Deconvolution" = "Deconv")
  d = mint[[regime]] %>% filter(method == m, pair == pn) %>%
    mutate(rna = reorder(rna, -rna_order), met = reorder(met, met_order), pair = short[[pn]])
  ggplot(d, aes(met, rna, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_point(data = filter(d, sig %in% TRUE, interaction > 0), aes(size = abs(interaction)), shape = 16, colour = "grey10") +
    geom_point(data = filter(d, sig %in% TRUE, interaction <= 0), aes(size = abs(interaction)), shape = 21, fill = "white", colour = "grey10", stroke = 0.35) +
    scale_size_area(limits = c(0, int_size_max), max_size = 2.6, oob = scales::squish, guide = "none") +
    scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                         limits = c(-lim_top[["row"]], lim_top[["row"]]), oob = scales::squish, guide = "none") +
    coord_cartesian(expand = FALSE) +
    labs(x = NULL, y = NULL) +
    ggh4x::facet_grid2(rows = vars(pair), switch = "y",
                       strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = unname(module_pal[[pair_family[[pn]]]]), colour = NA))) +
    ggh4x::force_panelsizes(rows = unit(nlevels(droplevels(d$rna)) * cell_de, "cm"), cols = unit(nlevels(droplevels(d$met)) * cell_de, "cm")) +
    theme_minimal(base_size = fs_lab) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = fs_lab),
          axis.text.y = if (first_col) element_text(size = fs_lab, hjust = 1) else element_blank(),
          axis.title = element_text(size = fs_lab, face = "italic"), axis.title.y = element_text(margin = margin(r = 1)),
          strip.placement = "outside",
          strip.text.y.left = if (first_col) element_text(angle = 90, face = "bold", size = fs_lab) else element_blank(),
          strip.background.y = if (first_col) NULL else element_blank(),
          plot.margin = margin(1, 1, 3, 1))
}
interaction_panel3 = function(regime, title) {
  methods = if (regime == "early") early_methods else late_methods
  pairs = names(pair_defs[[regime]])
  plots = unlist(lapply(pairs, function(pn) lapply(seq_along(methods), function(j) small_heat3(regime, methods[j], pn, j == 1))), recursive = FALSE)
  spacers = lapply(methods, function(m) ggplot() + theme_void() + ggh4x::force_panelsizes(cols = unit(0.1, "cm")))
  wrap_plots(c(lapply(methods, method_header), plots, spacers), ncol = length(methods), byrow = TRUE) +
    plot_annotation(title = title, subtitle = "rows: RNA option - columns: DNAm option",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = fs_title, margin = margin(b = 1)),
                                  plot.subtitle = element_text(hjust = 0.5, size = fs_lab, face = "italic", colour = "grey25", margin = margin(b = 4))))
}
int_early3 = interaction_panel3("early", "Early integration")
int_late3 = interaction_panel3("late", "Late integration")
# ---- legend block of B-E (drawn under D): colour scale, star of B-C, circles of D-E; coordinates in cm ----
pt_size = function(s) 2.6 * sqrt(pmin(s, int_size_max) / int_size_max) # same as scale_size_area(max_size = 2.6)
grad = data.frame(x = seq(0, 4, length.out = 81)) %>% mutate(v = (x / 4 * 2 - 1) * lim_top[["row"]])
ticks = data.frame(v = c(-0.1, -0.05, 0, 0.05, 0.1)) %>% mutate(x = (v / lim_top[["row"]] + 1) / 2 * 4)
sizes = data.frame(s = c(0.005, 0.01, 0.02, 0.03), lab = c("0.005", "0.01", "0.02", ">= 0.03")) %>% mutate(x = 4.2 + 0.95 * (row_number() - 1))
y1 = 2.35; y2 = 1.35; y3 = 0.55 # rows of the legend
p_leg = ggplot() +
  # B-E: colour scale
  annotate("text", x = 0, y = y1 + 0.42, label = "Score difference (B-E)",
           hjust = 0, size = fs_lab * pt, fontface = "bold") +
  geom_tile(data = grad, aes(x = x + 0.6, y = y1, fill = v), height = 0.22, width = 4 / 80) +
  geom_text(data = ticks, aes(x = x + 0.6, y = y1 - 0.25, label = sprintf("%+.2f", v)), vjust = 1, size = fs_lab * pt) +
  annotate("text", x = 0.6 + 4.25, y = y1, label = "(saturates at +/-0.10)", hjust = 0, size = fs_lab * pt, colour = "grey30") +
  # B, C: star
  annotate("text", x = 0, y = y2 + 0.22, label = "B, C:", hjust = 0, size = fs_lab * pt, fontface = "bold") +
  annotate("text", x = 0.75, y = y2 + 0.22 - 0.1, label = "*", size = 10 * pt, fontface = "bold", colour = "grey10") +
  annotate("text", x = 1.0, y = y2 + 0.22, label = "option different from the strategy median (Wilcoxon, BH-adjusted p < 0.05)", hjust = 0, size = fs_lab * pt) +
  # D, E: circles
  annotate("text", x = 0, y = y3 + 0.42, label = "D, E:", hjust = 0, size = fs_lab * pt, fontface = "bold") +
  annotate("text", x = 0.75, y = y3 + 0.42, label = "significant interaction contrast (Wilcoxon, BH-adjusted p < 0.05)", hjust = 0, size = fs_lab * pt) +
  annotate("point", x = 0.85, y = y3, shape = 16, size = 2, colour = "grey10") +
  annotate("text", x = 1.05, y = y3, label = "positive", hjust = 0, size = fs_lab * pt) +
  annotate("point", x = 2.15, y = y3, shape = 21, size = 2, fill = "white", colour = "grey10", stroke = 0.35) +
  annotate("text", x = 2.35, y = y3, label = "negative", hjust = 0, size = fs_lab * pt) +
  annotate("text", x = 0.75, y = y3 - 0.45, label = "size = abs(interaction contrast):", hjust = 0, size = fs_lab * pt) +
  geom_point(data = sizes, aes(x = x, y = y3 - 0.45, size = pt_size(s)), shape = 21, fill = "grey60", colour = "grey10", stroke = 0.3) +
  geom_text(data = sizes, aes(x = x + 0.1 + 0.05 * pt_size(s), y = y3 - 0.45, label = lab), hjust = 0, size = fs_lab * pt) +
  scale_size_identity() +
  scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                       limits = c(-lim_top[["row"]], lim_top[["row"]]), guide = "none") +
  coord_fixed(xlim = c(0, 8.6), ylim = c(-0.2, 3.0), expand = FALSE, clip = "off") +
  theme_void()
# ---- layout: A; B | C; (D over the legend) | E, D and E aligned at the top ----
top_align = function(p) wrap_plots(p, ggplot() + theme_void() + ggh4x::force_panelsizes(cols = unit(0.1, "cm")), ncol = 1) +
  plot_annotation(theme = theme(plot.background = element_blank())) # free-height spacer below; transparent: must not hide the axis title of A
row_BC = wrap_plots(wrap_elements(full = top_align(pB5)) + labs(tag = "B") + tagE, wrap_elements(full = top_align(pC5)) + labs(tag = "C") + tagE,
                    nrow = 1, widths = c(8.8, 9.2))
col_D = wrap_plots(wrap_elements(full = int_early3) + labs(tag = "D") + tagE, p_leg, ncol = 1, heights = c(6.4, 3.0))
row_DE = wrap_plots(col_D, wrap_elements(full = int_late3) + labs(tag = "E") + tagE, nrow = 1, widths = c(8.8, 9.2))

# ---- figure4_split -------------------------------------------------------------------------------------------------
# legend block (coordinates in cm): colour scale, then the star of the options or the circles of the interactions
legend_block = function(score_title, star = NULL, circles = NULL) {
  y1 = if (is.null(circles)) 0.95 else 2.35
  p = ggplot() +
    annotate("text", x = 0, y = y1 + 0.42, label = score_title, hjust = 0, size = fs_lab * pt, fontface = "bold") +
    geom_tile(data = grad, aes(x = x + 0.6, y = y1, fill = v), height = 0.22, width = 4 / 80) +
    geom_text(data = ticks, aes(x = x + 0.6, y = y1 - 0.25, label = sprintf("%+.2f", v)), vjust = 1, size = fs_lab * pt) +
    annotate("text", x = 0.6 + 4.25, y = y1, label = "(saturates at +/-0.10)", hjust = 0, size = fs_lab * pt, colour = "grey30")
  if (!is.null(star)) p = p +
    annotate("text", x = 9.6, y = y1 + 0.42, label = paste0(star, ":"), hjust = 0, size = fs_lab * pt, fontface = "bold") +
    annotate("text", x = 9.75, y = y1 - 0.1, label = "*", size = 10 * pt, fontface = "bold", colour = "grey10") +
    annotate("text", x = 10.0, y = y1, label = "option different from the strategy median\n(Wilcoxon, BH-adjusted p < 0.05)",
             hjust = 0, size = fs_lab * pt, lineheight = 0.9)
  if (!is.null(circles)) { y3 = 0.55; sizes = sizes %>% mutate(x = x + 0.35); p = p +
    annotate("text", x = 0, y = y3 + 0.42, label = paste0(circles, ":"), hjust = 0, size = fs_lab * pt, fontface = "bold") +
    annotate("text", x = 0.75, y = y3 + 0.42, label = "significant interaction contrast (Wilcoxon, BH-adjusted p < 0.05)", hjust = 0, size = fs_lab * pt) +
    annotate("point", x = 0.85, y = y3, shape = 16, size = 2, colour = "grey10") +
    annotate("text", x = 1.05, y = y3, label = "positive", hjust = 0, size = fs_lab * pt) +
    annotate("point", x = 2.15, y = y3, shape = 21, size = 2, fill = "white", colour = "grey10", stroke = 0.35) +
    annotate("text", x = 2.35, y = y3, label = "negative", hjust = 0, size = fs_lab * pt) +
    annotate("text", x = 0.75, y = y3 - 0.45, label = "size = abs(interaction contrast):", hjust = 0, size = fs_lab * pt) +
    geom_point(data = sizes, aes(x = x, y = y3 - 0.45, size = pt_size(s)), shape = 21, fill = "grey60", colour = "grey10", stroke = 0.3) +
    geom_text(data = sizes, aes(x = x + 0.1 + 0.05 * pt_size(s), y = y3 - 0.45, label = lab), hjust = 0, size = fs_lab * pt) }
  p + scale_size_identity() +
    scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                         limits = c(-lim_top[["row"]], lim_top[["row"]]), guide = "none") +
    coord_fixed(xlim = c(0, if (is.null(star)) 8.6 else 17.4), ylim = c(if (is.null(circles)) 0.45 else -0.2, y1 + 0.65), expand = FALSE, clip = "off") +
    theme_void()
}

# ---- supp_figure4_interactions -------------------------------------------------------------------------------------
# supplementary figure: the RNA x DNAm combinations (former D, E), aligned at the top; legend under A
col_A_supp = wrap_plots(wrap_elements(full = int_early3) + labs(tag = "A") + tagE,
                        legend_block("Score difference (A, B)", circles = "A, B"), ncol = 1, heights = c(6.4, 3.0))
supp_figure4_interactions = wrap_plots(col_A_supp, wrap_elements(full = int_late3) + labs(tag = "B") + tagE, nrow = 1, widths = c(8.8, 9.2))
ggsave(file.path(fig_dir, "SuppFig6.pdf"), plot = supp_figure4_interactions, width = 18, height = 10, units = "cm")
