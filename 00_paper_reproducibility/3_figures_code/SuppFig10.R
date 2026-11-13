# Supplementary Figure 10: best pipeline of each method re-optimised in each dataset.

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

# ---- supp_best_pipelines -------------------------------------------------------------------------------------------
supp_fontsize = 6 # A4 page: A above B, 18 cm wide
heat_A = best_pipeline_heatmap(ei_ds, early_methods, modules_early, lev_ei, main = "Early integration",
                               score_limits = heat_score_limits, hatch_keys = heat_hatch,
                               base_size = supp_fontsize, show = FALSE)
heat_B = best_pipeline_heatmap(li_ds, late_methods, modules_late, lev_li, extra_methods = uniomic_methods,
                               score_limits = heat_score_limits, hatch_keys = heat_hatch,
                               main = "Late integration and uni-omic pipelines",
                               base_size = supp_fontsize, show = FALSE)
supp_best = wrap_elements(full = heat_A) + labs(tag = "A") + tag_theme +
  wrap_elements(full = heat_B) + labs(tag = "B") + tag_theme + plot_layout(ncol = 1)
ggsave(file.path(fig_dir, "SuppFig10.pdf"), plot = supp_best, width = 18, height = 25, units = "cm")
