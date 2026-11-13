# Figure 6 of the article: generalist pipelines, benchmark and unseen data (lollipop version).

data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_04_figure6.R")

# ---- heatmaps_fun --------------------------------------------------------------------------------------------------
# one row per (integration method, dataset): its best pipeline. Columns: method and dataset (colour strips), score,
# then the module choices (black = level selected; hatched = level selected by the best pipeline of the dataset
# among `methods`). extra_methods: rows added below (e.g. the uni-omic pipelines), not used for the hatching; their
# modules that do not apply (nopp / nofs / node) are left blank. score_limits: limits of the colour scale of the
# score (default: range of the scores shown); hatch_keys: "method|dataset" of the rows to hatch (default: best of
# `methods` in each dataset). base_size: font size; file = NULL: no pdf; show = FALSE: not drawn. Returns the plot.
best_pipeline_heatmap = function(df, methods, modules, lev, file = NULL, main, base_size = 11, show = TRUE,
                                 width = 15, height = 12, extra_methods = NULL, score_limits = NULL, hatch_keys = NULL) {
  not_applicable = c("nopp", "nofs", "node")
  all_methods = c(methods, extra_methods)
  best = df %>%
    filter(method %in% all_methods) %>%
    group_by(method, dataset) %>%
    slice_max(score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(method = factor(method, levels = all_methods)) %>%
    arrange(method, dataset) %>%
    group_by(dataset) %>%
    mutate(best_of_dataset = if (is.null(hatch_keys)) method %in% methods & score == max(score[method %in% methods])
                             else paste(method, dataset, sep = "|") %in% hatch_keys) %>%
    ungroup() %>%
    # numeric row position inside the panel of its method (first dataset on top), used to draw the hatching
    group_by(method) %>% mutate(y = rev(seq_len(n()))) %>% ungroup()
  row_labels = paste0(method_label(as.character(best$method)), " | ", best$dataset)
  best$row = factor(row_labels, levels = rev(row_labels))

  # module choices: one column per level selected at least once (levels in the order of `lev`); the modules
  # that do not apply to a pipeline (nopp / nofs / node) give no column and are left blank on its row
  long = bind_rows(lapply(names(modules), function(mod) {
    vec = as.character(best[[modules[[mod]]]])
    vals = intersect(c(lev[[modules[[mod]]]], sort(unique(vec))), setdiff(unique(vec), not_applicable))
    expand.grid(i = seq_len(nrow(best)), x = seq_along(vals)) %>%
      mutate(module = mod, level = vals[x], selected = vec[i] == level, applies = !vec[i] %in% not_applicable)
  })) %>%
    mutate(row = best$row[i], y = best$y[i], method = best$method[i],
           hatched = selected & best$best_of_dataset[i],
           module = factor(module, levels = names(modules)),
           col = factor(paste(module, level, sep = "."), levels = unique(paste(module, level, sep = "."))))

  # hatching: diagonal segments inside each hatched cell (tile of width and height 1)
  h = 0.5
  hatch = long %>% filter(hatched) %>%
    tidyr::crossing(t = c(-0.6, -0.2, 0.2, 0.6)) %>%
    mutate(x0 = pmax(-h, -h - t), x1 = pmin(h, h - t),
           xs = x + x0, xe = x + x1, ys = y + x0 + t, ye = y + x1 + t)

  common_theme = theme_minimal(base_size = base_size) +
    theme(panel.grid = element_blank(), panel.spacing = unit(3, "pt"), strip.text.y = element_blank(),
          strip.text.x = element_text(face = "bold"), axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  # method strip (the dataset is written at the end of each row)
  annot = best %>% transmute(row, method, x = factor("method"), fill = unname(method_colors[as.character(method)]), top = "")
  p_annot = ggplot(annot, aes(x = x, y = row, fill = fill)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_identity() +
    facet_grid(method ~ top, scales = "free_y", space = "free_y") +
    common_theme + theme(axis.text.y = element_blank())

  # score: colour scale stretched over the range of the scores shown, so that the differences are visible
  p_score = ggplot(best %>% mutate(x = "", top = "Score"), aes(x = x, y = row, fill = score)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", score)), size = base_size * 0.3) +
    scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd")[1:7], limits = if (is.null(score_limits)) range(best$score) else score_limits,
                         guide = "none") +
    facet_grid(method ~ top, scales = "free_y", space = "free_y") +
    common_theme + theme(axis.text.y = element_blank())

  # module choices
  p_bin = ggplot(long, aes(x = col, y = row)) +
    geom_tile(data = filter(long, applies), aes(fill = selected & !hatched), colour = "grey60", linewidth = 0.25) +
    geom_tile(data = filter(long, !applies), fill = "white", colour = NA) + # module that does not apply: blank
    geom_segment(data = hatch, aes(x = xs, xend = xe, y = ys, yend = ye), linewidth = 0.35, inherit.aes = FALSE) +
    scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "black"), guide = "none") +
    scale_x_discrete(labels = function(x) sub("^[^.]*\\.", "", x)) +
    scale_y_discrete(position = "right") +
    facet_grid(method ~ module, scales = "free", space = "free",
               labeller = labeller(module = function(x) sub(" ", "\n", x))) + # two lines: long names were cut
    common_theme

  n_cols = nlevels(droplevels(long$col))
  p = wrap_plots(p_annot, p_score, p_bin, nrow = 1, widths = c(1, 3, n_cols)) +
    plot_annotation(title = main, theme = theme(plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5)))
  if (!is.null(file)) ggsave(file, plot = p, width = width, height = height)
  if (show) print(p)
  invisible(p)
}

# shared by all the heatmaps: colour scale of the score, and the row hatched in each dataset (its best pipeline
# over all the methods: early, late and uni-omic)
best_all = bind_rows(ei_ds %>% filter(method %in% early_methods),
                     li_ds %>% filter(method %in% c(late_methods, uniomic_methods))) %>%
  group_by(method, dataset) %>% summarise(score = max(score), .groups = "drop")
heat_score_limits = range(best_all$score)
heat_hatch = best_all %>% group_by(dataset) %>% slice_max(score, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(key = paste(method, dataset, sep = "|")) %>% pull(key)

# ---- recipes_fun ---------------------------------------------------------------------------------------------------
# the n_rec distinct generalist pipelines of method m with the best mean relative gain over the best uni-omic pipeline
gain_select = function(m) {
  core_cols = unname(unlist(rec_modules[[m]]))
  pipelines_long %>%
    filter(method == m, dataset %in% rank_datasets) %>%
    left_join(best_dna, by = "dataset") %>%
    left_join(best_rna, by = "dataset") %>%
    mutate(gain_rel = (score - pmax(best_onlyDNAm, best_onlyRNA)) / pmax(best_onlyDNAm, best_onlyRNA)) %>%
    group_by(id) %>%
    summarise(n_ds = n_distinct(dataset), mean_gain = mean(gain_rel), across(all_of(core_cols), first), .groups = "drop") %>%
    filter(n_ds == length(rank_datasets)) %>%
    arrange(desc(mean_gain)) %>%
    mutate(core = do.call(paste, c(lapply(across(all_of(core_cols)), as.character), sep = "|"))) %>%
    distinct(core, .keep_all = TRUE) %>%
    slice_head(n = n_rec) %>%
    mutate(top_rank = row_number(), row_label = as.character(top_rank))
}

# the four panels for the pipelines of `sel` (one row each: id, top_rank, row_label and the module columns)
# stats: summaries of panel 4 ("top": mean rank in %, "gain": mean gain over the best uni-omic pipeline)
rec_panels = function(m, tag, sel = rec_select(m), y_title = "recommended generalist pipeline (rank)", stats = c("top", "gain")) {
  mods = rec_modules[[m]]
  n_row = nrow(sel)

  # percentile rank (0 = best) of the retained pipelines in each dataset, and relative gain over the
  # best uni-omic pipeline of the dataset (the better of the best onlyDNAm and the best onlyRNA)
  d = pipelines_long %>%
    filter(method == m, dataset %in% rank_datasets) %>%
    group_by(dataset) %>%
    mutate(pct = rank(-score, ties.method = "average") / n()) %>%
    ungroup() %>%
    inner_join(select(sel, id, top_rank), by = "id", relationship = "many-to-many") %>% # a pipeline may fill several rows
    left_join(best_dna, by = "dataset") %>%
    left_join(best_rna, by = "dataset") %>%
    mutate(best_uniomic = pmax(best_onlyDNAm, best_onlyRNA),
           gain_rel = (score - best_uniomic) / best_uniomic,
           dataset = factor(dataset, levels = rank_datasets))

  # panel 1: the generalist pipeline
  recipe = bind_rows(lapply(names(mods), function(mod)
    data.frame(top_rank = sel$top_rank, module = mod, level = as.character(sel[[mods[[mod]]]]), stringsAsFactors = FALSE)))
  recipe$module = factor(recipe$module, levels = names(mods))
  p1 = ggplot(recipe, aes(x = module, y = top_rank, fill = level)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = level), size = 1.9) +
    scale_fill_manual(values = level_colors, na.value = "grey90", guide = "none") +
    scale_x_discrete(position = "top", labels = function(x) sub(" ", "\n", x)) + # two lines: long names overlapped
    scale_y_reverse(breaks = sel$top_rank, labels = sel$row_label, expand = c(0, 0)) +
    labs(title = bquote(bold(.(method_label(m)) * ": the best ") * bolditalic("generalist pipelines")), x = NULL, y = y_title, tag = tag) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6.5, face = "bold"),
          axis.text.y = element_text(size = 6.5), axis.title.y = element_text(size = 7),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5), plot.tag = element_text(face = "bold", size = 14))

  # panel 2: rank of the generalist pipeline in each dataset (top x % of the pipelines of its method)
  p2 = ggplot(d, aes(x = dataset, y = top_rank, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", 100 * pct)), size = 1.9) +
    scale_fill_gradientn(colours = rev(brewer.pal(7, "YlGnBu")), limits = c(0, 0.1), oob = scales::squish,
                         labels = scales::percent, name = "rank in the dataset (top x % of the pipelines of the method)",
                         guide = guide_colourbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"),
                                                 title.position = "top", title.hjust = 0.5)) +
    scale_x_discrete(position = "top") +
    scale_y_reverse(breaks = seq_len(n_row), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = " ") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 6.5, face = "bold"))

  # panel 3: relative gain over the best uni-omic pipeline of the dataset
  p3 = ggplot(d, aes(x = dataset, y = top_rank, fill = gain_rel)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.1f", 100 * gain_rel)), size = 1.9) +
    scale_fill_gradient2(low = div_cols[["low"]], mid = div_cols[["mid"]], high = div_cols[["high"]], midpoint = 0,
                         limits = c(-0.15, 0.15), oob = scales::squish, labels = scales::percent,
                         name = "relative gain over the best uni-omic pipeline of the dataset",
                         guide = guide_colourbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"),
                                                 title.position = "top", title.hjust = 0.5)) +
    scale_x_discrete(position = "top") +
    scale_y_reverse(breaks = seq_len(n_row), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = " ") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 6.5, face = "bold"))

  # panel 4: summaries
  summ = d %>% group_by(top_rank) %>%
    summarise(top = 100 * mean(pct), gain = 100 * mean(gain_rel), .groups = "drop") %>%
    tidyr::pivot_longer(-top_rank, names_to = "stat", values_to = "value") %>%
    filter(stat %in% stats) %>%
    mutate(stat = factor(stat, levels = c("top", "gain"), labels = c("mean\ntop %", "mean gain\nvs best\nuni-omic %")),
           label = ifelse(stat == "mean\ntop %", sprintf("%.1f", value), sprintf("%+.1f", value)))
  p4 = ggplot(summ, aes(x = stat, y = top_rank, label = label)) +
    geom_text(size = 2.1) +
    scale_x_discrete(position = "top") +
    # same vertical range as the tiles (no tile here: rows 1 and n_rec were cut); reversed scale: limits as c(max, min)
    scale_y_reverse(breaks = seq_len(n_row), limits = c(n_row + 0.5, 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = " ") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6.5, face = "bold", lineheight = 0.9))
  list(p1, p2, p3, p4)
}

# ---- panel_E -------------------------------------------------------------------------------------------------------
family_defs = list(
  "Cohort\ncorrelation" = c("pearson_tot_norm", "spearman_tot_norm", "pearson_col_norm", "spearman_col_norm"),
  "Cell-type\ncorrelation" = c("pearson_row_norm", "spearman_row_norm"),
  "Error\nmetrics" = c("rmse_norm", "mae_norm", "aitchison_norm", "jsd_norm"),
  "Geometric\nmetrics" = c("aid_norm", "sdid_norm"))
top_methods = fig6_methods
# generalist pipelines re-run with the function blocks (05_rerun_top_pipelines.R): those of figure 6, and those of ConcatStd
# and Mean (used by 12_exploration_other_methods.Rmd, 13_celltype_rna_vs_dnam.Rmd, 15_hrRMSE_fixed_recipes.Rmd and the check of
# 00_data_processing.Rmd)
rerun_methods = c(top_methods, "concatscale", "limean")
top1_all = bind_rows(lapply(rerun_methods, function(m) gain_select(m) %>% slice_head(n = 1) %>% transmute(method = m, id)))
top1 = top1_all %>% filter(method %in% top_methods)
# modules of the generalist pipelines, for their re-run with the function blocks (05_rerun_top_pipelines.R)
top1_modules = pipelines_long %>% inner_join(top1_all, by = c("method", "id")) %>% distinct(method, id, .keep_all = TRUE) %>%
  select(method, any_of(c(all_functions_ei, all_functions_li)), id)
# block names of the framework in the file (read by the re-run scripts 05_*.R): "SumNorm" and "wRLR" written back as "Scale" and "RLRpoisson"
write.csv(top1_modules %>% mutate(across(starts_with("preprocessing"), ~ ifelse(.x == "SumNorm", "Scale", as.character(.x))),
                                 across(starts_with("deconvolution"), ~ ifelse(.x == "wRLR", "RLRpoisson", as.character(.x))),
                                 id = gsub("wRLR", "RLRpoisson", gsub("SumNorm", "Scale", id))),
          file.path(out_dir, "top1_pipelines.csv"), row.names = FALSE)
perf = pipelines_long %>% inner_join(top1, by = c("method", "id"))
# geometric mean of the metrics available (JSD is missing in VITR and SDN4: their error family uses RMSE, MAE, Aitchison)
geo_mean = function(d, cols) { g = exp(rowMeans(log(as.matrix(d[, cols])), na.rm = TRUE)); g[is.nan(g)] = NA; g }
fam_long = bind_rows(lapply(names(family_defs), function(f) {
  data.frame(method = perf$method, dataset = perf$dataset, score = perf$score, family = f,
             value = geo_mean(perf, family_defs[[f]]))
})) %>%
  mutate(method = factor(method, levels = top_methods), family = factor(family, levels = names(family_defs)))

# check: the geometric mean of the 4 families gives back the stored aggregated score
complete_ds = perf %>% filter(if_all(all_of(unlist(family_defs)), ~ !is.na(.x))) %>% distinct(method, dataset)
chk = fam_long %>% semi_join(complete_ds, by = c("method", "dataset")) %>%
  group_by(method, dataset) %>% summarise(score = first(score), rebuilt = exp(mean(log(value))), .groups = "drop")

top_label = setNames(paste0(method_label(top_methods), "\n(top 1)"), top_methods)
top_cols = method_colors[top_methods] # the method colours of all the figures (points with a dark outline)
# aggregated score (5th panel): the stored score_aggreg of the pipeline in the dataset
# (VIVO is not shown in the aggregated score panel: its score uses a different formula)
fam_all = bind_rows(fam_long, fam_long %>% distinct(method, dataset, score) %>% filter(dataset != "VIVO") %>% mutate(family = "Aggregated\nscore", value = score)) %>%
  mutate(family = factor(as.character(family), levels = c(names(family_defs), "Aggregated\nscore")))
fam_plot = fam_all %>% filter(!is.na(value)) %>% mutate(x = as.integer(method), is_vivo = dataset == "VIVO")
fam_med = fam_plot %>% group_by(family, method, x) %>% summarise(med = median(value), .groups = "drop")

# one panel per family; one point per dataset, the points of a dataset linked across the methods (grey line); bar = median
p_E = ggplot(fam_plot, aes(x = x, y = value)) +
  geom_line(aes(group = dataset), colour = "grey85", linewidth = 0.3) +
  geom_segment(data = fam_med, aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med), colour = "black", linewidth = 0.7) +
  geom_point(aes(fill = method, shape = is_vivo), size = 2, stroke = 0.3, colour = "grey15") +
  facet_wrap(~ family, nrow = 1) +
  scale_fill_manual(values = top_cols, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = c(`FALSE` = "dataset", `TRUE` = "VIVO (different score)"), name = NULL) +
  scale_x_continuous(breaks = seq_along(top_methods), labels = top_label, expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  labs(title = "Top pipelines: score by family of metrics (one point per dataset, bar = median)", x = NULL, y = "Family score") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 6.5, lineheight = 0.9, angle = 45, hjust = 1, vjust = 1), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

# ---- panel_F -------------------------------------------------------------------------------------------------------
rerun_name = c(OnlyMet = "onlyDNAm", OnlyRna = "onlyRNA", concatscale = "concatscale", concatnoscale = "concatnoscale",
               limean = "limean", limeanRMSE = "limeanRMSE", TunedJ = "TunedJ")
rerun_file = file.path(data_dir, "reruns_top1_recipes_scores.csv") # 05_results/<pipeline>/scores.csv of the five generalist pipelines, gathered by 05c_export_reruns.R
p_F = NULL
if (file.exists(rerun_file)) {
  rr = read.csv(rerun_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(method %in% top_methods, mix == 2) %>%
    mutate(dataset = factor(dataset_levels[match(sub("_mix2$", "", dataset), dataset_levels)], levels = dataset_levels))
  fam_F = bind_rows(lapply(names(family_defs), function(f)
    data.frame(method = rr$method, dataset = rr$dataset, score = rr$score_aggreg, family = f, value = geo_mean(rr, family_defs[[f]]))))
  fam_F = bind_rows(fam_F, fam_F %>% distinct(method, dataset, score) %>% filter(dataset != "VIVO") %>% mutate(family = "Aggregated\nscore", value = score)) %>% # VIVO not in the aggregated score panel
    mutate(method = factor(method, levels = top_methods), family = factor(family, levels = c(names(family_defs), "Aggregated\nscore")))
  write.csv(fam_F %>% mutate(family = gsub("\n", " ", family)), file.path(out_dir, "supp_figure6_metric_families_unseen.csv"), row.names = FALSE)
  d_F = fam_F %>% filter(!is.na(value)) %>% mutate(x = as.integer(method), is_vivo = dataset == "VIVO")
  med_F = d_F %>% group_by(family, method, x) %>% summarise(med = median(value), .groups = "drop")
  p_F = ggplot(d_F, aes(x = x, y = value)) +
    geom_line(aes(group = dataset), colour = "grey85", linewidth = 0.3) +
    geom_segment(data = med_F, aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med), colour = "black", linewidth = 0.7) +
    geom_point(aes(fill = method, shape = is_vivo), size = 2, stroke = 0.3, colour = "grey15") +
    facet_wrap(~ family, nrow = 1) +
    scale_fill_manual(values = top_cols, guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = c(`FALSE` = "dataset", `TRUE` = "VIVO (different score)"), name = NULL) +
    scale_x_continuous(breaks = seq_along(top_methods), labels = top_label, expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1), oob = scales::squish) +
    labs(title = "Top pipelines run on unseen datasets: score by family of metrics (one point per dataset, bar = median)",
         x = NULL, y = "Family score") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 6.5, lineheight = 0.9, angle = 45, hjust = 1, vjust = 1), strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  p_F
} else cat("re-run scores not found (run 05_rerun_top_pipelines.R): panel F not drawn\n")

# ---- figure7_tests -------------------------------------------------------------------------------------------------
agg_scores = bind_rows(
  fam_all %>% filter(family == "Aggregated\nscore") %>% transmute(panel = "A", method = as.character(method), dataset = as.character(dataset), score = value),
  if (!is.null(p_F)) fam_F %>% filter(family == "Aggregated\nscore") %>% transmute(panel = "B", method = as.character(method), dataset = as.character(dataset), score = value))
paired_tests = function(d) {
  ref = d %>% filter(method == "OnlyMet") %>% select(dataset, ref = score)
  # only the strategies shown in panels B and D (onlyRNA is not shown, so not tested): BH over the 3 integration strategies
  d %>% filter(method %in% c(fig6_early, fig6_late, "TunedJ")) %>% inner_join(ref, by = "dataset") %>%
    group_by(method) %>%
    summarise(n = n(), n_better = sum(score > ref), median_diff = median(score - ref),
              p = wilcox.test(score, ref, paired = TRUE, exact = TRUE, alternative = "greater")$p.value, .groups = "drop") %>%
    mutate(p_BH = p.adjust(p, method = "BH"))
}
tests_fig7 = agg_scores %>% group_by(panel) %>% group_modify(~ paired_tests(.x)) %>% ungroup() %>%
  mutate(method = factor(method, levels = top_methods)) %>% arrange(panel, method)

# ---- scheme_fun ----------------------------------------------------------------------------------------------------
dataset_files = c(invitro = "VITR", insilicopseudobulk = "SBN5", insilicodirichletNoDep = "SDN5",
                  insilicodirichletNoDep4CTsource = "SDN4", insilicodirichletNoDep6CTsource = "SDN6",
                  insilicodirichletEMFA = "SDE5", insilicodirichletEMFAImmuneLowProp = "SDEL",
                  insilicodirichletCopule = "SDC5", invivo = "VIVO") # file name of the mixes -> short name
# schematic of the two top pipelines: one row per omic branch, boxes = modules, arrows = data flow
t1 = rename_pp(read.csv(file.path(out_dir, "top1_pipelines.csv"), stringsAsFactors = FALSE)) # "SumNorm" in the figures
tj = t1[t1$method == "TunedJ", ]; om = t1[t1$method == "OnlyMet", ]
box = function(x, y, label, fill, w = 1.9) data.frame(x = x, y = y, label = label, fill = fill, w = w)
step_cols = c(data = "grey90", pp = unname(module_pal[["preprocessing"]]), fs = unname(module_pal[["feature_selection"]]),
              de = unname(module_pal[["deconvolution"]]), li = unname(module_pal[["late_integration"]]), out = "white")
branch = function(y, omic, r, side) {
  pp = if (side == "rna") r$preprocessing_RNA else r$preprocessing_MET
  fs = if (side == "rna") r$feature_selection_RNA else r$feature_selection_MET
  de = if (side == "rna") r$deconvolution_rna else r$deconvolution_met
  bind_rows(box(1, y, paste0(omic, "\nmixture + reference"), step_cols[["data"]]),
            box(3.1, y, paste0("Pre-processing\n", pp), step_cols[["pp"]]),
            box(5.2, y, paste0("Feature selection\n", fs), step_cols[["fs"]]),
            box(7.3, y, paste0("Deconvolution\n", de), step_cols[["de"]]))
}
scheme = function(r, title, rna) {
  b = bind_rows(if (rna) branch(2, "RNA", r, "rna"), branch(1, "DNAm", r, "met"))
  ymid = if (rna) 1.5 else 1
  li = if (rna) box(10.3, ymid, "Late integration (PriorMean)\nclassic, basal: DNAm only\nendo, fibro, immune:\nmean of RNA and DNAm", step_cols[["li"]], w = 3.3) else
    box(10.3, ymid, "onlyDNAm\n(DNAm prediction)", step_cols[["li"]], w = 3.3)
  out = box(12.9, ymid, "cell type\nproportions", step_cols[["out"]], w = 1.3)
  arrows = bind_rows(
    b %>% group_by(y) %>% arrange(x) %>% transmute(x0 = x + w / 2, x1 = lead(x) - lead(w) / 2, y0 = y, y1 = y) %>% filter(!is.na(x1)) %>% ungroup(),
    b %>% filter(x == 7.3) %>% transmute(x0 = x + w / 2, x1 = li$x - li$w / 2, y0 = y, y1 = ymid),
    data.frame(x0 = li$x + li$w / 2, x1 = out$x - out$w / 2, y0 = ymid, y1 = ymid))
  all = bind_rows(b, li, out)
  ggplot() +
    geom_segment(data = arrows, aes(x = x0, xend = x1, y = y0, yend = y1), arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.35, colour = "grey30") +
    geom_tile(data = all, aes(x = x, y = y, width = w, height = ifelse(grepl("PriorMean", label), 0.95, 0.72), fill = fill),
              colour = "grey40", linewidth = 0.3) +
    geom_text(data = all, aes(x = x, y = y, label = label), size = 2, lineheight = 0.85) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(0.55, 2.45)) + scale_x_continuous(limits = c(0, 13.6)) +
    labs(title = title) + theme_void(base_size = 8) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
stars = function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
fam_levels = c(names(family_defs), "Aggregated\nscore")
brackets = function(panel_id) {
  tests_fig7 %>% filter(panel == panel_id) %>% arrange(method) %>%
    transmute(family = factor("Aggregated\nscore", levels = fam_levels), x0 = 1, x1 = as.integer(method),
              y = 0.99 + 0.045 * (row_number() - 1), label = stars(p_BH), signif = p_BH < 0.05)
}
# (y axis extended so that the brackets sit above the points)
add_p = function(p, panel_id) {
  b = brackets(panel_id)
  suppressMessages(p +
    geom_segment(data = b, aes(x = x0, xend = x1, y = y, yend = y), linewidth = 0.3, inherit.aes = FALSE) +
    geom_segment(data = b, aes(x = x0, xend = x0, y = y, yend = y - 0.012), linewidth = 0.3, inherit.aes = FALSE) +
    geom_segment(data = b, aes(x = x1, xend = x1, y = y, yend = y - 0.012), linewidth = 0.3, inherit.aes = FALSE) +
    geom_text(data = b, aes(x = (x0 + x1) / 2, y = y + ifelse(signif, 0.002, 0.012), label = label),
              size = ifelse(b$signif, 3.2, 2), vjust = 0, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0.4, max(b$y) + 0.05), breaks = seq(0.4, 1, 0.1), oob = scales::squish))
}
cap = labs(caption = "One-sided paired Wilcoxon signed-rank test, method > onlyDNAm (8 datasets, VIVO excluded), Benjamini-Hochberg adjusted: * p < 0.05, ** p < 0.01, *** p < 0.001, ns: not significant")
# number of features of the reference after the preprocessing and after the feature selection of each top 1 generalist pipeline
# (counted on the files of the re-runs by 05c_export_reruns.R)
n_feat_tab = if (file.exists(file.path(data_dir, "reruns_top1_recipes_n_features.csv"))) read.csv(file.path(data_dir, "reruns_top1_recipes_n_features.csv"), stringsAsFactors = FALSE)
n_feat = function(omic, pp, fs = NULL) {
  pp = ifelse(pp == "SumNorm", "Scale", pp) # the table keeps the name of the block
  if (is.null(n_feat_tab)) return(NA)
  r = n_feat_tab[n_feat_tab$omic == omic & n_feat_tab$preprocessing == pp, ]
  if (!is.null(fs)) r = r[r$feature_selection == fs, ]
  if (nrow(r) == 0) return(NA)
  format(if (is.null(fs)) r$n_after_preprocessing[1] else r$n_after_feature_selection[1], big.mark = ",")
}
branch_n = function(y, omic, r, side) { # branch() with the numbers of features in the preprocessing and selection boxes
  pp = if (side == "rna") r$preprocessing_RNA else r$preprocessing_MET
  fs = if (side == "rna") r$feature_selection_RNA else r$feature_selection_MET
  om = if (side == "rna") "RNA" else "MET"
  bl = branch(y, omic, r, side)
  bl$label[bl$x == 3.1] = paste0("Pre-processing\n", pp, "\n", n_feat(om, pp), if (side == "rna") " genes" else " CpGs")
  bl$label[bl$x == 5.2] = paste0("Feature selection\n", fs, "\n", n_feat(om, pp, fs), if (side == "rna") " genes" else " CpGs")
  bl
}
# schematic of any top pipeline: uni-omic, late integration (two arms) or early integration (concatenation, then one deconvolution)
scheme_any = function(m) {
  branch = branch_n
  r = t1[t1$method == m, ]
  late_txt = c(limean = "Late integration (Mean)\nmean of RNA and DNAm",
               limeanRMSE = "Late integration (wMean)\nweighted mean of RNA and DNAm\n(weights: inverse\nreconstruction errors)",
               TunedJ = "Late integration (PriorMean)\nclassic, basal: DNAm only\nendo, fibro, immune:\nmean of RNA and DNAm")
  if (m %in% c("OnlyMet", "OnlyRna")) {
    side = if (m == "OnlyMet") "met" else "rna"
    b = branch(1.5, if (side == "met") "DNAm" else "RNA", r, side); ymid = 1.5 # one arm: centred vertically
    mid = box(10.3, ymid, paste0(method_label(m), "\n(", if (side == "met") "DNAm" else "RNA", " prediction)"), step_cols[["li"]], w = 3.3)
    arrows_last = b %>% filter(x == 7.3) %>% transmute(x0 = x + w / 2, x1 = mid$x - mid$w / 2, y0 = y, y1 = ymid)
    extra = NULL
  } else if (m %in% names(late_txt)) {
    b = bind_rows(branch(2, "RNA", r, "rna"), branch(1, "DNAm", r, "met")); ymid = 1.5
    mid = box(10.3, ymid, late_txt[[m]], step_cols[["li"]], w = 3.3)
    arrows_last = b %>% filter(x == 7.3) %>% transmute(x0 = x + w / 2, x1 = mid$x - mid$w / 2, y0 = y, y1 = ymid)
    extra = NULL
  } else { # early integration: the two pre-processed and selected omics are concatenated, then deconvolved once
    b = bind_rows(branch(2, "RNA", r, "rna"), branch(1, "DNAm", r, "met")) %>% filter(x < 7)
    ymid = 1.5
    ei = box(7.6, ymid, paste0("Early integration\n(", m, ")"), unname(module_pal[["early_integration"]]), w = 2.1)
    mid = box(10.3, ymid, paste0("Deconvolution\n", r$deconvolution), step_cols[["de"]], w = 2.4)
    arrows_last = bind_rows(b %>% filter(x == 5.2) %>% transmute(x0 = x + w / 2, x1 = ei$x - ei$w / 2, y0 = y, y1 = ymid),
                            data.frame(x0 = ei$x + ei$w / 2, x1 = mid$x - mid$w / 2, y0 = ymid, y1 = ymid))
    extra = ei
  }
  out = box(12.9, ymid, "cell type\nproportions", step_cols[["out"]], w = 1.3)
  arrows = bind_rows(
    b %>% group_by(y) %>% arrange(x) %>% transmute(x0 = x + w / 2, x1 = lead(x) - lead(w) / 2, y0 = y, y1 = y) %>% filter(!is.na(x1)) %>% ungroup(),
    arrows_last, data.frame(x0 = mid$x + mid$w / 2, x1 = out$x - out$w / 2, y0 = ymid, y1 = ymid))
  all = bind_rows(b, extra, mid, out)
  ggplot() +
    geom_segment(data = arrows, aes(x = x0, xend = x1, y = y0, yend = y1), arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.35, colour = "grey30") +
    geom_tile(data = all, aes(x = x, y = y, width = w, height = ifelse(grepl("PriorMean", label), 0.96, 0.84), fill = fill),
              colour = "grey40", linewidth = 0.3) +
    geom_text(data = all, aes(x = x, y = y, label = label), size = 1.9, lineheight = 0.8) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) + # compact: no margin around the boxes + scale_x_continuous(limits = c(0, 13.6)) +
    labs(title = bquote(bold(.(method_label(m)) * ": top 1 of the best ") * bolditalic("generalist pipelines"))) +
    theme_void(base_size = 8) + theme(plot.title = element_text(hjust = 0.5, colour = unname(top_cols[m]), margin = margin(b = 1)),
                                      plot.margin = margin(1, 4, 1, 4))
}
p_scheme_A = wrap_plots(scheme_any("OnlyMet") + labs(tag = "A") + tag_theme, scheme_any("OnlyRna"),
                        scheme_any(fig6_early), scheme_any(fig6_late),
                        scheme_any("TunedJ"), plot_spacer(), ncol = 2)

# ---- panels_BD -----------------------------------------------------------------------------------------------------
fixed_rows = c(fig6_early, fig6_late, "TunedJ") # rows of B and D (first on top); onlyRNA left out: it adds little
cond_lab = c(A = "Benchmark datasets", B = "Unseen datasets")
# aggregated score of the top 1 generalist pipelines in the 9 datasets (VIVO included: the relative gain is computed within a dataset)
scores_9 = bind_rows(fam_long %>% distinct(method, dataset, score) %>% mutate(panel = "A"),
                     fam_F %>% distinct(method, dataset, score) %>% mutate(panel = "B")) %>%
  mutate(method = as.character(method), dataset = as.character(dataset))
gain_fixed = scores_9 %>%
  group_by(panel, dataset) %>%
  mutate(ref = score[method == "OnlyMet"], gain = (score - ref) / ref) %>%
  ungroup() %>%
  filter(method %in% fixed_rows)
gain_fixed = bind_rows(gain_fixed, gain_fixed %>% filter(dataset != "VIVO") %>% group_by(panel, method) %>%
                         summarise(gain = mean(gain), .groups = "drop") %>% mutate(dataset = "mean")) %>%
  mutate(dataset = factor(dataset, levels = c(dataset_levels, "mean")),
         block = factor(ifelse(dataset == "mean", "mean", "dataset"), levels = c("dataset", "mean")),
         group = factor(ifelse(method == "OnlyRna", "reference", "integration"), levels = c("reference", "integration")), # small gap: onlyRNA is not an integration
         method = factor(method, levels = fixed_rows)) # facet rows: first level on top

lim_gain = 0.06 # colour limit (squished beyond: the values are written in the cells)
gain_cols = div_cols # same diverging palette as the other heatmaps
gain_heat = function(panel_id) {
  prime = if (panel_id == "B") "'" else ""
  ggplot(filter(gain_fixed, panel == panel_id), aes(x = dataset, y = method, fill = gain)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(abs(gain) < 0.0005, "0.0", sprintf("%+.1f", 100 * gain)), fontface = ifelse(block == "mean", "bold", "plain"),
                  colour = abs(gain) > 0.75 * lim_gain), size = 2.2) +
    geom_vline(data = data.frame(block = factor("dataset", levels = c("dataset", "mean"))), aes(xintercept = length(rank_datasets) + 0.5),
               linetype = "dashed", colour = "grey30", linewidth = 0.4) + # VIVO apart: different score
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    scale_fill_gradient2(low = gain_cols[["low"]], mid = gain_cols[["mid"]], high = gain_cols[["high"]], midpoint = 0,
                         limits = c(-lim_gain, lim_gain), oob = scales::squish, labels = function(x) sprintf("%+.0f %%", 100 * x),
                         name = "relative gain over the onlyDNAm generalist pipeline (red: better than DNAm; blue: worse; colours saturate at +/-6 %)",
                         guide = guide_colourbar(barwidth = unit(6, "cm"), barheight = unit(0.25, "cm"), title.position = "top", title.hjust = 0.5)) +
    scale_x_discrete(labels = function(x) ifelse(x == "mean", "mean", paste0(x, ifelse(x == "VIVO", "*", ""), prime))) +
    # rows: one strip per strategy, of the colour of the method (as in figure 4)
    ggh4x::facet_grid2(method ~ block, scales = "free", space = "free", switch = "y", labeller = labeller(method = function(x) method_label(x)),
                       strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = unname(method_colors[fixed_rows]), colour = NA),
                                                   text_y = ggh4x::elem_list_text(colour = text_on(unname(method_colors[fixed_rows]))))) +
    labs(title = cond_lab[[panel_id]], x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), strip.text.x = element_blank(), panel.spacing.x = unit(4, "pt"), panel.spacing.y = unit(1, "pt"),
          strip.placement = "outside", strip.text.y.left = element_text(size = 7, face = "bold", angle = 0, hjust = 1, margin = margin(1, 3, 1, 3)),
          axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9))
}
p_B = gain_heat("A")
p_C = gain_heat("B")

# ---- panels_CE -----------------------------------------------------------------------------------------------------
fam_levels = c(names(family_defs), "Aggregated\nscore")
fam_methods = setdiff(top_methods, "OnlyRna") # strategies of panels C and E (onlyRNA left out: it adds little)
family_tests = function(fam, cond) {
  d = fam %>% filter(dataset != "VIVO", !is.na(value)) %>% mutate(method = as.character(method), family = as.character(family))
  ref = d %>% filter(method == "OnlyMet") %>% select(family, dataset, ref = value)
  d %>% filter(method %in% fam_methods, method != "OnlyMet") %>% inner_join(ref, by = c("family", "dataset")) %>% # onlyRNA not tested: not shown
    group_by(family, method) %>%
    summarise(n = n(), n_better = sum(value > ref), median_diff = median(value - ref),
              p = wilcox.test(value, ref, paired = TRUE, exact = TRUE, alternative = "greater")$p.value, .groups = "drop") %>%
    group_by(family) %>% mutate(p_BH = p.adjust(p, method = "BH")) %>% ungroup() %>%
    mutate(condition = cond)
}
tests_fam = bind_rows(family_tests(fam_all, "benchmark datasets"), if (!is.null(p_F)) family_tests(fam_F, "unseen datasets"))

# dotted line at the median of onlyDNAm, and the stars of the tests above each method, in every family
add_fam_tests = function(p, cond, med) {
  ref_line = med %>% filter(method == "OnlyMet") %>% select(family, med) %>%
    mutate(family = factor(as.character(family), levels = fam_levels))
  st = tests_fam %>% filter(condition == cond) %>%
    transmute(family = factor(family, levels = fam_levels), x = match(method, top_methods), label = stars(p_BH))
  st = bind_rows(st, data.frame(family = factor(fam_levels, levels = fam_levels), x = 1, label = "ref."))
  suppressMessages(p +
    geom_hline(data = ref_line, aes(yintercept = med), linetype = "dotted", colour = "grey20", linewidth = 0.45) +
    geom_text(data = st, aes(x = x, y = 1.015, label = label), size = ifelse(st$label %in% c("ns", "ref."), 1.9, 3), vjust = 0, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0.4, 1.06), breaks = seq(0.4, 1, 0.1), oob = scales::squish) +
    scale_x_continuous(breaks = seq_along(top_methods), labels = method_label(top_methods), expand = expansion(add = 0.5)) + # no "(top 1)": labels overlapped
    theme(axis.text.x = element_text(size = 6.5, angle = 50, hjust = 1, vjust = 1)))
}
p_Cfam = add_fam_tests(p_E, "benchmark datasets", fam_med) + labs(title = "Benchmark datasets: score by family of metrics")
p_Efam = if (!is.null(p_F)) add_fam_tests(p_F, "unseen datasets", med_F) + labs(title = "Unseen datasets: score by family of metrics")

# ---- figure6 -------------------------------------------------------------------------------------------------------
pt = 1 / ggplot2::.pt # geom_text size for 1 pt (geom_text sizes are in mm)
fs_min = 6            # smallest font size (pt)
tag6 = theme(plot.tag = element_text(face = "bold", size = 14))

# ---- A: the generalist pipelines, one line per omic arm, columns = modules ----
fmt_n = function(omic, pp, fs = NULL) {
  n = n_feat(if (omic == "RNA") "RNA" else "MET", pp, fs)
  paste(n, if (omic == "RNA") "genes" else "CpGs")
}
recipe_rows = bind_rows(lapply(top_methods, function(m) {
  r = t1[t1$method == m, ]
  arms = switch(m, OnlyMet = "DNAm", OnlyRna = "RNA", c("RNA", "DNAm"))
  bind_rows(lapply(arms, function(o) {
    pp = if (o == "RNA") r$preprocessing_RNA else r$preprocessing_MET
    fsel = if (o == "RNA") r$feature_selection_RNA else r$feature_selection_MET
    de = if (m == fig6_early) NA else paste0("Deconvolution\n", if (o == "RNA") r$deconvolution_rna else r$deconvolution_met)
    data.frame(method = m, omic = o, pp = paste0(pp, "\n", fmt_n(o, pp)), fs = paste0(fsel, "\n", fmt_n(o, pp, fsel)), de = de)
  }))
}))
integ_txt = c(OnlyMet = "No integration\n(DNAm prediction)", OnlyRna = "No integration\n(RNA prediction)",
              setNames(paste0("Deconvolution: ", t1$deconvolution[t1$method == fig6_early], "\n(one joint deconvolution of the concatenated data)"), fig6_early),
              limean = "Late integration\nmean of the RNA and DNAm predictions",
              limeanRMSE = "Late integration\nmean of the RNA and DNAm predictions weighted by\nthe inverse of their reconstruction errors",
              TunedJ = "Late integration\nclassic, basal: DNAm only;\nendo, fibro, immune: mean of RNA and DNAm")
# vertical layout: one line per arm, a gap between generalist pipelines (first generalist pipeline on top)
gap = 0.35
recipe_rows = recipe_rows %>% mutate(i = row_number(), block = match(method, top_methods),
                                     y = -(i + gap * (block - 1)))
blocks = recipe_rows %>% group_by(method, block) %>% summarise(ymid = mean(y), ymin = min(y), ymax = max(y), .groups = "drop")
# columns (x centre, width)
cx = c(name = 0.95, omic = 2.55, pp = 4.35, fs = 6.85, de = 9.25, int = 13.15)
cw = c(omic = 1.0, pp = 2.25, fs = 2.4, de = 1.9, int = 5.3)
h = 0.94 # tile height (one arm)
tiles = bind_rows(
  recipe_rows %>% transmute(x = cx[["omic"]], w = cw[["omic"]], y, hh = h, label = omic, fill = "grey90"),
  recipe_rows %>% transmute(x = cx[["pp"]], w = cw[["pp"]], y, hh = h, label = pp, fill = step_cols[["pp"]]),
  recipe_rows %>% transmute(x = cx[["fs"]], w = cw[["fs"]], y, hh = h, label = fs, fill = step_cols[["fs"]]),
  recipe_rows %>% filter(!is.na(de)) %>% transmute(x = cx[["de"]], w = cw[["de"]], y, hh = h, label = de, fill = step_cols[["de"]]),
  # early integration (Concat): the concatenation comes before the deconvolution: integration in the deconvolution
  # column, joint deconvolution of the concatenated data in the last column (both over the two arms)
  blocks %>% filter(method == fig6_early) %>%
    transmute(x = cx[["de"]], w = cw[["de"]], y = ymid, hh = ymax - ymin + h,
              label = c(concatscale = "Early integration\n(scaled\nconcatenation)", concatnoscale = "Early integration\n(concatenation)")[[fig6_early]],
              fill = unname(module_pal[["early_integration"]])),
  # last column: late integration (olive), joint deconvolution of the early integration (copper), nothing to integrate for the uni-omic pipelines (white)
  blocks %>% transmute(x = cx[["int"]], w = cw[["int"]], y = ymid, hh = ymax - ymin + h, label = integ_txt[method],
                       fill = ifelse(method == fig6_early, step_cols[["de"]], ifelse(method %in% uniomic_methods, "white", step_cols[["li"]]))))
arrows_A = bind_rows(
  recipe_rows %>% transmute(x0 = cx[["omic"]] + cw[["omic"]] / 2, x1 = cx[["pp"]] - cw[["pp"]] / 2, y0 = y, y1 = y),
  recipe_rows %>% transmute(x0 = cx[["pp"]] + cw[["pp"]] / 2, x1 = cx[["fs"]] - cw[["fs"]] / 2, y0 = y, y1 = y),
  recipe_rows %>% transmute(x0 = cx[["fs"]] + cw[["fs"]] / 2, x1 = cx[["de"]] - cw[["de"]] / 2, y0 = y, y1 = y),
  recipe_rows %>% left_join(blocks, by = c("method", "block")) %>% filter(method != fig6_early) %>%
    transmute(x0 = cx[["de"]] + cw[["de"]] / 2, x1 = cx[["int"]] - cw[["int"]] / 2, y0 = y, y1 = ymid),
  blocks %>% filter(method == fig6_early) %>% transmute(x0 = cx[["de"]] + cw[["de"]] / 2, x1 = cx[["int"]] - cw[["int"]] / 2, y0 = ymid, y1 = ymid))
# headers of the two first columns only: the deconvolution and the integration are named in their boxes (their order differs for the early integration)
heads = data.frame(x = cx[c("pp", "fs")], label = c("Pre-processing", "Feature selection"))
p_A6 = ggplot() +
  geom_segment(data = arrows_A, aes(x = x0, xend = x1, y = y0, yend = y1), arrow = arrow(length = unit(1, "mm"), type = "closed"),
               linewidth = 0.3, colour = "grey35") +
  geom_tile(data = tiles, aes(x = x, y = y, width = w, height = hh, fill = fill), colour = "grey45", linewidth = 0.25) +
  geom_text(data = tiles, aes(x = x, y = y, label = label), size = fs_min * pt, lineheight = 0.8,
            colour = text_on(tiles$fill)) + # text colour from the box fill
  geom_text(data = heads, aes(x = x, y = 0.05, label = label), size = 6.5 * pt, fontface = "bold", vjust = 0) +
  # strategy names on bands of the colour of the method, as in figure 4
  geom_tile(data = blocks, aes(x = cx[["name"]], y = ymid, height = ymax - ymin + h, fill = unname(method_colors[method])), width = 1.75, colour = NA) +
  geom_text(data = blocks, aes(x = cx[["name"]], y = ymid, label = method_label(method)), size = 6.5 * pt, fontface = "bold",
            colour = text_on(unname(method_colors[blocks$method]))) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, cx[["int"]] + cw[["int"]] / 2 + 0.05), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.45))) +
  labs(tag = "A", title = "Best generalist pipeline of each strategy: one combination of modules, applied unchanged to every dataset") +
  theme_void() + tag6 + theme(plot.margin = margin(2, 2, 2, 2), plot.title = element_text(face = "bold", hjust = 0.5, size = 7.5, margin = margin(b = 3)))

# ---- B, D: relative gain over the onlyDNAm generalist pipeline ----
heat6 = function(panel_id) {
  prime = if (panel_id == "B") "'" else "" # ASCII: the unicode prime is not rendered in this locale
  d = filter(gain_fixed, panel == panel_id) %>%
    mutate(lab = ifelse(abs(gain) < 0.0005, "0.0", ifelse(abs(gain) >= 0.1, sprintf("%+.0f", 100 * gain), sprintf("%+.1f", 100 * gain))))
  ggplot(d, aes(x = dataset, y = method, fill = gain)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab, fontface = ifelse(block == "mean", "bold", "plain"), colour = abs(gain) > 0.75 * lim_gain), size = fs_min * pt) +
    geom_vline(data = data.frame(block = factor("dataset", levels = c("dataset", "mean"))), aes(xintercept = length(rank_datasets) + 0.5),
               linetype = "dashed", colour = "grey30", linewidth = 0.35) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
    scale_fill_gradient2(low = gain_cols[["low"]], mid = gain_cols[["mid"]], high = gain_cols[["high"]], midpoint = 0,
                         limits = c(-lim_gain, lim_gain), oob = scales::squish, labels = function(x) sprintf("%+.0f %%", 100 * x),
                         name = "B, D: relative gain over the onlyDNAm generalist pipeline\n(colours saturate at +/-6 %)",
                         guide = guide_colourbar(barwidth = unit(3.5, "cm"), barheight = unit(0.2, "cm"), title.position = "top", title.hjust = 0.5)) +
    scale_x_discrete(labels = function(x) ifelse(x == "mean", "mean", paste0(x, ifelse(x == "VIVO", "*", ""), prime))) +
    # rows: one strip per strategy, of the colour of the method (as in figure 4)
    ggh4x::facet_grid2(method ~ block, scales = "free", space = "free", switch = "y", labeller = labeller(method = function(x) method_label(x)),
                       strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = unname(method_colors[fixed_rows]), colour = NA),
                                                   text_y = ggh4x::elem_list_text(colour = text_on(unname(method_colors[fixed_rows]))))) +
    labs(title = cond_lab[[panel_id]], x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(), strip.text.x = element_blank(), panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(1, "pt"),
          strip.placement = "outside", strip.text.y.left = element_text(size = 6.5, face = "bold", angle = 0, hjust = 1, margin = margin(1, 3, 1, 3)),
          axis.text.x = element_text(size = fs_min, angle = 50, hjust = 1, vjust = 1), axis.text.y = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 7), legend.title = element_text(size = fs_min), legend.text = element_text(size = fs_min))
}

# ---- C, E: score by family of metrics, tests against onlyDNAm ----
fam6 = function(d, cond) {
  d = d %>% filter(method %in% fam_methods) %>%
    mutate(family = factor(as.character(family), levels = fam_levels), x = match(as.character(method), fam_methods))
  med = d %>% group_by(family, method, x) %>% summarise(med = median(value), .groups = "drop")
  ref_line = med %>% filter(method == "OnlyMet")
  st = tests_fam %>% filter(condition == cond) %>%
    transmute(family = factor(family, levels = fam_levels), x = match(method, fam_methods), label = stars(p_BH))
  st = bind_rows(st, data.frame(family = factor(fam_levels, levels = fam_levels), x = 1, label = "ref."))
  ggplot(d, aes(x = x, y = value)) +
    geom_hline(data = ref_line, aes(yintercept = med), linetype = "dotted", colour = "grey20", linewidth = 0.4) +
    geom_line(aes(group = dataset), colour = "grey82", linewidth = 0.25) +
    geom_segment(data = med, aes(x = x - 0.32, xend = x + 0.32, y = med, yend = med), colour = "black", linewidth = 0.55) +
    geom_point(aes(fill = method, shape = is_vivo), size = 1.1, stroke = 0.2, colour = "grey15") +
    geom_text(data = st, aes(x = x, y = 1.02, label = label), size = ifelse(st$label %in% c("ns", "ref."), fs_min * pt, 9 * pt), vjust = 0, inherit.aes = FALSE) +
    facet_wrap(~ family, nrow = 1) +
    scale_fill_manual(values = top_cols, guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = c(`FALSE` = "C, E: dataset", `TRUE` = "VIVO (different score)"), name = NULL) +
    scale_x_continuous(breaks = seq_along(fam_methods), labels = method_label(fam_methods), expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0.4, 1.08), breaks = seq(0.4, 1, 0.1), oob = scales::squish) +
    labs(title = paste0(cond_lab[[c(`benchmark datasets` = "A", `unseen datasets` = "B")[[cond]]]], ": score by family of metrics"), x = NULL, y = "Family score") +
    theme_bw(base_size = 7) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.spacing.x = unit(2, "pt"),
          axis.text.x = element_text(size = fs_min, angle = 55, hjust = 1, vjust = 1), axis.text.y = element_text(size = fs_min),
          strip.text = element_text(face = "bold", size = fs_min, lineheight = 0.9, margin = margin(1.5, 1, 1.5, 1)),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 7), legend.text = element_text(size = fs_min))
}

# B and D in their own frames: their rotated labels and strips do not shift the axes of C and E
row_BC6 = wrap_plots(wrap_elements(full = heat6("A") + guides(fill = "none")) + labs(tag = "B") + tag6,
                     fam6(fam_plot, "benchmark datasets") + labs(tag = "C") + tag6 + guides(shape = "none"),
                     nrow = 1, widths = c(1, 1.35))
row_DE6 = wrap_plots(wrap_elements(full = heat6("B") + theme(legend.position = "bottom", legend.box.spacing = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0))) +
                       labs(tag = "D") + tag6,
                     fam6(d_F, "unseen datasets") + labs(tag = "E") + tag6,
                     nrow = 1, widths = c(1, 1.35), guides = "collect") &
  theme(legend.position = "bottom", legend.box.spacing = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0))
# heights in cm: A keeps the height of its boxes, the rows B-C and D-E are compact
figure6 = wrap_plots(p_A6, row_BC6, row_DE6, ncol = 1, heights = c(8.9, 5.6, 6.5))

# ---- figure6_lollipop ----------------------------------------------------------------------------------------------
lim_lolli = 0.10 # axis cut (relative gain)
lolli6 = function(panel_id) {
  prime = if (panel_id == "B") "'" else ""
  d = filter(gain_fixed, panel == panel_id) %>%
    mutate(dataset = factor(as.character(dataset), levels = rev(c(dataset_levels, "mean"))),
           g = pmax(pmin(gain, lim_lolli), -lim_lolli), cut = abs(gain) > lim_lolli,
           sign = factor(ifelse(gain >= 0, "better", "worse"), levels = c("better", "worse")))
  ggplot(d, aes(y = dataset)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_hline(data = data.frame(block = factor("dataset", levels = c("dataset", "mean"))), aes(yintercept = 1.5),
               linetype = "dotted", colour = "grey30", linewidth = 0.35) + # VIVO apart: different score
    geom_segment(aes(x = 0, xend = g, yend = dataset, colour = sign), linewidth = 0.6) +
    geom_point(aes(x = g, fill = sign, shape = block), size = ifelse(d$block == "mean", 1.9, 1.4), stroke = 0.25, colour = "grey15") +
    geom_text(data = filter(d, cut), aes(x = g, label = sprintf("%+.0f", 100 * gain), hjust = ifelse(gain > 0, 1.35, -0.35)),
              size = fs_min * pt, vjust = -0.45) +
    scale_colour_manual(values = c(better = div_cols[["high"]], worse = div_cols[["low"]]), guide = "none") +
    scale_fill_manual(values = c(better = div_cols[["high"]], worse = div_cols[["low"]]), guide = "none") +
    scale_shape_manual(values = c(dataset = 21, mean = 23), guide = "none") +
    scale_x_continuous(limits = c(-lim_lolli, lim_lolli) * 1.05, breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
                       labels = function(x) ifelse(abs(x) > 0.09, "", ifelse(x == 0, "0", sprintf("%+.0f", 100 * x))), expand = c(0, 0)) + # +/-10: grid only (labels of adjacent facets overlap)
    scale_y_discrete(labels = function(x) ifelse(x == "mean", "mean", paste0(x, ifelse(x == "VIVO", "*", ""), prime))) +
    ggh4x::facet_grid2(block ~ method, scales = "free_y", space = "free_y", labeller = labeller(method = function(x) method_label(x)),
                       strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = unname(method_colors[fixed_rows]), colour = NA),
                                                   text_x = ggh4x::elem_list_text(colour = text_on(unname(method_colors[fixed_rows])), face = rep("bold", length(fixed_rows))))) +
    labs(title = cond_lab[[panel_id]], x = "relative gain over the onlyDNAm generalist pipeline (%)", y = NULL) +
    theme_bw(base_size = 7) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(2, "pt"),
          strip.text.y = element_blank(), strip.text.x = element_text(size = 6.5, margin = margin(1.5, 1, 1.5, 1)),
          axis.text = element_text(size = fs_min), axis.title.x = element_text(size = fs_min),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 7))
}
row_BC6l = wrap_plots(lolli6("A") + labs(tag = "B") + tag6,
                      fam6(fam_plot, "benchmark datasets") + labs(tag = "C") + tag6 + guides(shape = "none"),
                      nrow = 1, widths = c(1, 1.35))
row_DE6l = wrap_plots(lolli6("B") + labs(tag = "D") + tag6,
                      fam6(d_F, "unseen datasets") + labs(tag = "E") + tag6,
                      nrow = 1, widths = c(1, 1.35)) &
  theme(legend.position = "bottom", legend.box.spacing = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0))
figure6_lollipop = wrap_plots(p_A6, row_BC6l, row_DE6l, ncol = 1, heights = c(8.9, 5.8, 6.3))
ggsave(file.path(fig_dir, "Fig6.pdf"), plot = figure6_lollipop, width = 18, height = 21, units = "cm", device = cairo_pdf)
