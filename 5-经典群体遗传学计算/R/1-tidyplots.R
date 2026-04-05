library(tidyverse)
library(tidyplots)

setwd("/mnt/l/26-2026年4月3日临时项目/5-77SNP/script/")

# output_filter_dir <- "../output/5-过滤"
# dir.create(output_filter_dir, recursive = TRUE, showWarnings = FALSE)

df_within <- read_tsv('../output/1-within/within_group_summary.tsv')


p1 <- df_within |>
  tidyplot(
    x = group,
    y = n_samples,
    color = group
  ) |>
  add_mean_bar() |>
  adjust_x_axis(rotate_labels = 90) |>
  sort_x_axis_labels(.reverse = TRUE) |>
  add_data_labels(label = n_samples, label_position = "above") |>
  adjust_title("Number of samples") 


df_boots_within <- read_tsv("../output/4-bootstrap/within_group_bootstrap_replicates.tsv")

df_boots_within |>select(metric) |> unique()

# df_boots_within_filtered <- bind_rows(
#   df_boots_within |>
#     filter(metric == "pi", value <= 0.005),
#   df_boots_within |>
#     filter(metric == "theta_w", value <= 0.01),
#   df_boots_within |>
#     filter(metric == "tajima_d", value >= -2.2),
#   df_boots_within |>
#     filter(metric == "haplotype_diversity")
# )

# write_tsv(
#   df_boots_within_filtered,
#   file.path(output_filter_dir, "within_group_bootstrap_replicates.filtered.tsv")
# )

# df_boots_within_filtered_stats <- df_boots_within_filtered |>
#   group_by(metric, group) |>
#   summarise(
#     n = n(),
#     mean = mean(value, na.rm = TRUE),
#     sd = sd(value, na.rm = TRUE),
#     min = min(value, na.rm = TRUE),
#     q1 = quantile(value, 0.25, na.rm = TRUE),
#     median = median(value, na.rm = TRUE),
#     q3 = quantile(value, 0.75, na.rm = TRUE),
#     max = max(value, na.rm = TRUE),
#     .groups = "drop"
#   )

write_tsv(
  df_boots_within_filtered_stats,
  file.path(output_filter_dir, "within_group_bootstrap_replicates.filtered.summary_stats.tsv")
)

p4 <- df_boots_within |>
  filter(metric == 'pi') |>
  tidyplot(
    x = group,
    y = value,
    color = group
  ) |>
  add_data_points_beeswarm(size = 0.5,  alpha = 0.05, rasterize = TRUE, rasterize_dpi = 600) |>
  add_violin() |>
  adjust_x_axis(rotate_labels = 90) |>
  remove_legend() |>
  sort_x_axis_labels(.reverse = TRUE)  |>
  adjust_title("Bootstrap of nucleotide diversity (pi)")

p5 <- df_boots_within |>
  filter(metric == 'theta_w') |>
  tidyplot(
    x = group,
    y = value,
    color = group
  ) |>
  add_data_points_beeswarm(size = 0.5,  alpha = 0.05, rasterize = TRUE, rasterize_dpi = 600) |>
  add_violin() |>
  adjust_x_axis(rotate_labels = 90) |>
  remove_legend() |>
  sort_x_axis_labels(.reverse = TRUE)  |>
  adjust_title("Bootstrap of Watterson's theta (theta_w)")

p6 <- df_boots_within |>
  filter(metric == 'tajima_d') |>
  tidyplot(
    x = group,
    y = value,
    color = group
  ) |>
  add_data_points_beeswarm(size = 0.5,  alpha = 0.05,rasterize = TRUE, rasterize_dpi = 600) |>
  add_violin() |>
  adjust_x_axis(rotate_labels = 90) |>
  remove_legend() |>
  sort_x_axis_labels(.reverse = TRUE) |>
  adjust_title("Bootstrap of Tajima's D")

p7 <- df_boots_within |>
  filter(metric == 'haplotype_diversity') |>
  tidyplot(
    x = group,
    y = value,
    color = group
  ) |>
  add_data_points_beeswarm(size = 0.5,  alpha = 0.05, rasterize = TRUE, rasterize_dpi = 600) |>
  add_violin() |>
  adjust_x_axis(rotate_labels = 90) |>
  remove_legend() |>
  sort_x_axis_labels(.reverse = TRUE) |>
  adjust_title("Bootstrap of haplotype diversity")


#! 使用patchwork v1.3.0组合tidyplot对象需要使用如下代码
p_all <- patchwork::wrap_plots(p1, p4, p5, p6, p7) + 
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(filename = "../output/merge.pdf", plot = p_all, width = 20, height = 40, units = "cm")
