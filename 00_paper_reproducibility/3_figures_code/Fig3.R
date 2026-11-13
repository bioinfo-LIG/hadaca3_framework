# Figure 3 of the article: performance of the integration methods.


data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_01_figure3.R")

# ---- panel_A -------------------------------------------------------------------------------------------------------
half_violin = function(data, x_var, title, frame_col = "grey60", score = "score_aggreg", level_order,
                       show_y = TRUE, base_size = 9, label_size = 3, adjust = 1.3, frame_width = 1) {
  data[[x_var]] = droplevels(data[[x_var]])
  lv = level_order[level_order %in% levels(data[[x_var]])]
  data[[x_var]] = factor(as.character(data[[x_var]]), levels = lv)

  qfun = function(l, p) quantile(data[[score]][data[[x_var]] == l], probs = p, names = FALSE, na.rm = TRUE)
  sm = data.frame(level = factor(lv, levels = lv), i = seq_along(lv),
                  p05 = sapply(lv, qfun, p = 0.05), q1 = sapply(lv, qfun, p = 0.25), med = sapply(lv, qfun, p = 0.5),
                  q3 = sapply(lv, qfun, p = 0.75), p95 = sapply(lv, qfun, p = 0.95))
  rownames(sm) = NULL

  hv = do.call(rbind, lapply(seq_along(lv), function(i) {
    x = data[[score]][data[[x_var]] == lv[i]]
    x = x[!is.na(x)]
    d = density(x, adjust = adjust, from = max(min(x), 0), to = min(max(x), 1), n = 256)
    dens = d$y / max(d$y)
    data.frame(level = lv[i], y = c(d$x, rev(d$x)), x = c(rep(i + 0.05, length(d$x)), rev(i + 0.05 + 0.42 * dens)))
  }))
  hv$level = factor(hv$level, levels = lv)

  ggplot() +
    geom_polygon(data = hv, aes(x = x, y = y, group = level, fill = level), colour = "grey30", linewidth = 0.25) +
    geom_segment(data = sm, aes(x = i - 0.12, xend = i - 0.12, y = p05, yend = p95), colour = "grey30", linewidth = 0.4) +
    geom_segment(data = sm, aes(x = i - 0.12, xend = i - 0.12, y = q1, yend = q3),
                 colour = "grey10", linewidth = 2.2 * base_size / 12, lineend = "butt") +
    geom_point(data = sm, aes(x = i - 0.12, y = med), shape = 21, fill = "white", colour = "grey10",
               size = 2.6 * base_size / 12, stroke = 0.6) +
    geom_text(data = sm, aes(x = i, y = 1.04, label = sprintf("%.2f", med)), fontface = "bold", size = label_size, vjust = 0) +
    scale_fill_manual(values = method_colors[lv]) +
    scale_x_continuous(breaks = seq_along(lv), labels = method_label(lv), limits = c(0.4, length(lv) + 0.6)) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    labs(title = title, x = NULL, y = "Aggregated score") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.4),
      axis.line = element_line(colour = "grey40", linewidth = 0.4),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size),
      panel.border = element_rect(colour = frame_col, fill = NA, linewidth = frame_width),
      legend.position = "none"
    ) +
    (if (show_y) NULL else theme(axis.title.y = element_blank(), axis.text.y = element_blank()))
}

violin_strip = function(res, base_size = 9, label_size = 3) {
  pA = half_violin(res$ei, "early_integration", "Early integration", "grey60",
                   level_order = early_methods, show_y = TRUE, base_size = base_size, label_size = label_size) +
    labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0.02, 0.98))
  pB = half_violin(droplevels(res$li[res$li$late_integration %in% late_methods, ]), "late_integration",
                   "Late integration", "grey60",
                   level_order = late_methods, show_y = FALSE, base_size = base_size, label_size = label_size)
  pC = half_violin(droplevels(res$li[res$li$late_integration %in% uniomic_methods, ]), "late_integration",
                   "Uni-omic", "grey60",
                   level_order = uniomic_methods, show_y = FALSE, base_size = base_size, label_size = label_size)
  wrap_plots(pA, pB, pC, nrow = 1, widths = c(5.2, 3.2, 2.2))
}

p_A = violin_strip(res_with)

# ---- panel_B -------------------------------------------------------------------------------------------------------
median_all = all_ds %>%
  group_by(dataset, method) %>%
  summarise(median_score = median(score), .groups = "drop")

# 95% bootstrap confidence interval of the median (2000 resamplings of the pipelines); saved, as it takes a few minutes
{
  set.seed(1)
  boot_ci = function(x, n_boot = 2000) quantile(replicate(n_boot, median(sample(x, replace = TRUE))), c(0.025, 0.975), names = FALSE)
  ci_all = all_ds %>%
    group_by(dataset, method) %>%
    summarise(ci = list(boot_ci(score)), .groups = "drop") %>%
    mutate(ci_lower = sapply(ci, `[`, 1), ci_upper = sapply(ci, `[`, 2)) %>%
    select(-ci)
}

median_all = median_all %>%
  left_join(ci_all, by = c("dataset", "method")) %>%
  mutate(method = factor(method, levels = method_order), is_vivo = dataset == "VIVO") %>%
  group_by(dataset) %>%
  mutate(label = ifelse(median_score == max(median_score), sprintf("%.2f", median_score), "")) %>%
  ungroup()

p_B = ggplot(median_all, aes(x = dataset, y = median_score, fill = method)) +
  geom_vline(xintercept = length(dataset_levels) - 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_col(aes(alpha = is_vivo), position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, group = method), position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 0.3, colour = "black") +
  geom_text(aes(y = ci_upper, label = label, group = method), position = position_dodge(width = 0.9),
            angle = 90, hjust = -0.15, vjust = 0.4, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.55), guide = "none") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.3)) +
  labs(x = NULL, y = "Median aggregated score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 1)))

# ---- panel_C -------------------------------------------------------------------------------------------------------
met_cols = c("preprocessing_mixMET", "feature_selection_mixMET", "preprocessing_MET", "feature_selection_MET")
# exact two-sided Wilcoxon signed-rank p-value, by enumeration of the 2^n sign permutations (n <= 15 here): zero
# differences are dropped and ties get mid-ranks (wilcox.test(exact = TRUE) falls back to the normal approximation
# as soon as there are ties or zeros)
exact_signrank_p = function(x) {
  x = x[x != 0]; n = length(x)
  if (n == 0) return(1)
  r = rank(abs(x)); t_obs = sum(r[x > 0]); e = sum(r) / 2
  signs = as.matrix(expand.grid(rep(list(c(0, 1)), n)))
  t_all = drop(signs %*% r)
  mean(abs(t_all - e) >= abs(t_obs - e) - 1e-9)
}
summarise_config = function(d) {
  d %>%
    mutate(config_key = do.call(paste, c(lapply(across(all_of(met_cols)), as.character), sep = "||"))) %>%
    group_by(dataset, config_key, method) %>%
    summarise(score = median(score), .groups = "drop")
}
onlymet_summary = summarise_config(all_ds %>% filter(method == "OnlyMet"))
int_summary = summarise_config(all_ds %>% filter(method %in% c(early_methods, late_methods)))

paired_vs_onlymet = int_summary %>%
  inner_join(onlymet_summary %>% select(dataset, config_key, score_ref = score), by = c("dataset", "config_key")) %>%
  group_by(dataset, method) %>%
  group_modify(function(d, key) {
    if (nrow(d) < 5) return(data.frame(n_pairs = nrow(d), estimate = NA, ci_lower = NA, ci_upper = NA, p_value = NA))
    wt = suppressWarnings(wilcox.test(d$score, d$score_ref, paired = TRUE, exact = FALSE, conf.int = TRUE, conf.level = 0.95))
    data.frame(n_pairs = nrow(d), estimate = unname(wt$estimate), ci_lower = wt$conf.int[1], ci_upper = wt$conf.int[2],
               p_value = exact_signrank_p(d$score - d$score_ref))
  }) %>%
  group_by(dataset) %>%
  mutate(p_adj_bonf = p.adjust(p_value, method = "bonferroni")) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      is.na(p_adj_bonf) | p_adj_bonf >= 0.05 ~ "n.s.",
      estimate < 0 & p_adj_bonf < 0.01 ~ "INT < onlyDNAm, p < 0.01",
      estimate < 0 ~ "INT < onlyDNAm, p < 0.05",
      p_adj_bonf < 0.01 ~ "INT > onlyDNAm, p < 0.01",
      TRUE ~ "INT > onlyDNAm, p < 0.05"),
    sig = factor(sig, levels = c("INT < onlyDNAm, p < 0.01", "INT < onlyDNAm, p < 0.05", "n.s.",
                                 "INT > onlyDNAm, p < 0.05", "INT > onlyDNAm, p < 0.01")),
    method = factor(method, levels = c(early_methods, late_methods)),
    # same order of the datasets as in the other figures (VIVO last)
    dataset = factor(as.character(dataset), levels = dataset_levels))

p_C = ggplot(paired_vs_onlymet, aes(x = estimate, y = method, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 0.4) +
  geom_point(size = 1.6) +
  facet_wrap(~ dataset, ncol = 3) +
  scale_colour_manual(values = c("INT < onlyDNAm, p < 0.01" = "#762A83", "INT < onlyDNAm, p < 0.05" = "#C2A5CF", "n.s." = "grey65",
                                 "INT > onlyDNAm, p < 0.05" = "#A6DBA0", "INT > onlyDNAm, p < 0.01" = "#1B7837"),
                      drop = FALSE, name = "Significance\n(Bonferroni-corrected)") +
  scale_y_discrete(labels = method_label) +
  scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.2)) + # few breaks: the labels of neighbouring facets overlapped
  labs(title = "Median paired difference vs. onlyDNAm (integration - DNAm)", x = NULL, y = "Integration strategies (INT)") +
  theme_bw() +
  theme(strip.text = element_text(size = 8), axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0, size = 10.5), legend.position = "bottom") + # left: the centred title ran into the tag of panel D
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# ---- panel_D -------------------------------------------------------------------------------------------------------
rank_median = median_all %>%
  group_by(dataset) %>%
  mutate(rank = rank(-median_score, ties.method = "average")) %>%
  ungroup() %>%
  mutate(x = as.integer(dataset))

p_D = ggplot(rank_median, aes(x = x, y = rank, colour = method, group = method)) +
  geom_vline(xintercept = length(dataset_levels) - 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_line(data = filter(rank_median, x < length(dataset_levels)), linewidth = 0.6) +
  geom_line(data = filter(rank_median, x >= length(dataset_levels) - 1), linewidth = 0.5, linetype = "dotted") +
  geom_point(aes(shape = dataset == "VIVO"), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_colour_manual(values = method_colors, breaks = method_order, labels = method_label, guide = "none") +
  scale_x_continuous(breaks = seq_along(dataset_levels), labels = dataset_levels) +
  scale_y_reverse(breaks = seq_along(method_order)) +
  labs(x = NULL, y = "Median rank (1 = best)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# ---- figure3 -------------------------------------------------------------------------------------------------------
pt = 1 / ggplot2::.pt
# one style for the four panels: thin grey frames, horizontal grid lines, short bold titles, tags at the same place
style3 = theme(panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.3), axis.line = element_blank(),
               panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
               strip.background = element_rect(fill = "grey92", colour = "grey50", linewidth = 0.3), strip.text = element_text(size = 6.5, face = "bold", margin = margin(1.5, 0, 1.5, 0)),
               axis.text = element_text(size = 6), axis.title = element_text(size = 6.5),
               plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5),
               plot.tag = element_text(face = "bold", size = 14), plot.tag.position = "topleft", # own space: the tags do not overlap the panels
               legend.text = element_text(size = 6), legend.title = element_text(size = 6.5), legend.key.size = unit(0.28, "cm"),
               legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))
# significance of C: red = better than onlyDNAm, blue = worse, as in the heatmaps of the other figures
sig_cols = c("INT < onlyDNAm, p < 0.01" = "#762A83", "INT < onlyDNAm, p < 0.05" = "#C2A5CF", "n.s." = "grey65",
             "INT > onlyDNAm, p < 0.05" = "#A6DBA0", "INT > onlyDNAm, p < 0.01" = "#1B7837")

# A: score distribution by method (three groups: early, late, uni-omic): full violins with a boxplot inside, median written above
dA = bind_rows(res_with$ei %>% transmute(method = as.character(early_integration), score = score_aggreg),
               res_with$li %>% transmute(method = as.character(late_integration), score = score_aggreg)) %>%
  filter(method %in% method_order, !is.na(score)) %>%
  mutate(method = factor(method, levels = method_order),
         group = factor(ifelse(method %in% early_methods, "Early integration", ifelse(method %in% late_methods, "Late integration", "Uni-omic")),
                        levels = c("Early integration", "Late integration", "Uni-omic")))
medA = dA %>% group_by(method, group) %>% summarise(m = median(score), .groups = "drop")
pA3 = ggplot(dA, aes(method, score, fill = method)) +
  geom_violin(scale = "width", width = 0.85, linewidth = 0.25, colour = "grey25") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", colour = "grey15", linewidth = 0.3) +
  geom_text(data = medA, aes(x = method, y = 1.03, label = sprintf("%.2f", m)), size = 6 * pt, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_discrete(labels = method_label) +
  scale_y_continuous(limits = c(0, 1.07), breaks = seq(0, 1, 0.25)) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  labs(title = "Score of the pipelines by method", tag = "A", x = NULL, y = "Aggregated score") +
  theme_bw(base_size = 7) + style3 + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5))

# B: median score by dataset (95 % bootstrap interval)
pB3 = suppressMessages(p_B + theme_bw(base_size = 7) +
  geom_text(aes(y = ci_upper, label = label, group = method), position = position_dodge(width = 0.9),
            angle = 90, hjust = -0.15, vjust = 0.4, size = 6 * pt, fontface = "bold") +
  scale_fill_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL) +
  labs(title = "Median score by dataset", tag = "B") + style3 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 1))))
pB3$layers[[4]] = NULL # the label layer of p_B (size 2.6 mm), replaced by the one above (6 pt)

# C: paired difference with onlyDNAm
pC3 = suppressMessages(p_C + theme_bw(base_size = 7) +
  # legend order: p < 0.05 before p < 0.01 in both rows (INT < and INT >)
  scale_colour_manual(values = sig_cols, drop = FALSE, name = "Significance\n(Bonferroni-corrected)",
                      breaks = c("INT < onlyDNAm, p < 0.05", "INT < onlyDNAm, p < 0.01", "n.s.",
                                 "INT > onlyDNAm, p < 0.05", "INT > onlyDNAm, p < 0.01")) +
  labs(title = "Paired difference vs. onlyDNAm", tag = "C") + style3 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)))

# D: rank of the median score
# points with a dark outline: the palest method colours stay visible
n_ds = length(dataset_levels)
pD3 = ggplot(rank_median, aes(x = x, y = rank, colour = method, group = method)) +
  geom_vline(xintercept = n_ds - 0.5, linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_line(data = filter(rank_median, x < n_ds), linewidth = 0.55) +
  geom_line(data = filter(rank_median, x >= n_ds - 1), linewidth = 0.45, linetype = "dotted") +
  geom_point(aes(fill = method, shape = dataset == "VIVO"), colour = "grey30", size = 1.4, stroke = 0.3) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
  scale_colour_manual(values = method_colors, guide = "none") + scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_continuous(breaks = seq_along(dataset_levels), labels = dataset_levels) +
  scale_y_reverse(breaks = seq_along(method_order)) +
  labs(title = "Rank of the median score", tag = "D", x = NULL, y = "Median rank (1 = best)") +
  theme_bw(base_size = 7) + style3 + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5))

figure3 = wrap_plots(pA3, pB3, wrap_plots(pC3, pD3, nrow = 1, widths = c(1.25, 1)),
                     ncol = 1, heights = c(6.4, 5.4, 9.2)) # heights in cm
ggsave(file.path(fig_dir, "Fig3.pdf"), plot = figure3, width = 18, height = 22, units = "cm")
