# Supplementary Figure 9: leaderboard-decay curves (top 100 pipelines of each method and dataset).

data_dir = "../00_results" # tables of the benchmark (see ../00_results/README.md)
fig_dir = "figures"
dir.create(fig_dir, showWarnings = FALSE)
source("setup_03_figure5.R")

# ---- supp_leaderboard_data -----------------------------------------------------------------------------------------
leaderboard = bind_rows(ei_ds, li_ds) %>%
  select(dataset, method, score) %>%
  group_by(dataset, method) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  mutate(rank = row_number(), n_total = n(), rank_pct = 100 * rank / n_total) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = method_order), dataset = factor(as.character(dataset), levels = dataset_levels))

leaderboard_plot = function(d, x_var, x_lab, title) {
  ggplot(d, aes(x = .data[[x_var]], y = score, colour = method)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ dataset, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = method_colors, breaks = method_order, labels = method_label, name = NULL, guide = "none") +
    labs(x = x_lab, y = "Aggregated score", title = title) +
    theme_bw(base_size = 7) +
    theme(legend.position = "bottom", strip.text = element_text(size = 7), axis.text = element_text(size = 6.5),
          panel.grid.minor = element_blank(), legend.key.height = unit(8, "pt"),
          plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5)) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(linewidth = 1)))
}

# ---- supp_leaderboard_top100 ---------------------------------------------------------------------------------------
supp_leaderboard_top100 = leaderboard_plot(filter(leaderboard, rank <= 100), "rank", "Rank among the evaluated pipelines (top 100)",
                                           "Leaderboard decay (top 100)")
ggsave(file.path(fig_dir, "SuppFig9.pdf"), plot = supp_leaderboard_top100, width = 18, height = 15, units = "cm") # A4 page: 18 cm wide
