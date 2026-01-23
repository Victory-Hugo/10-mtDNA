library(tidyverse)
library(tidyplots)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: Rscript vis.R <输入目录> <输出PDF路径>")
}

input_dir <- args[1]
out_file <- args[2]

paths <- list(
  dxy_csv = file.path(input_dir, "pairwise", "dxy.csv"),
  fst_pair = file.path(input_dir, "pairwise", "fst.csv"),
  pi_csv = file.path(input_dir, "population", "pi.csv"),
  tajima_d_csv = file.path(input_dir, "population", "tajima_d.csv"),
  theta_w_csv = file.path(input_dir, "population", "theta_w.csv"),
  population_count_csv = file.path(input_dir, "population", "population_counts.csv")
)

my_color <- c("#d55e00", "#f5c710", "#009e73", "#56b4e9", "#0072b2")

format_heatmap <- function(df, value_col, value_label = NULL) {
  df |>
    tidyplot(
      x = population1,
      y = population_2,
      fill = {{ value_col }}
    ) |>
    add_heatmap(rasterize = FALSE) |>
    adjust_colors(
      new_colors = new_color_scheme(my_color),
      breaks = scales::pretty_breaks(n = 5),
      labels = if (is.null(value_label)) {
        waiver()
      } else {
        value_label
      }
    )
}

format_bar <- function(df, y_col, label_vec) {
  df |>
    tidyplot(
      x = population,
      y = {{ y_col }},
      color = {{ y_col }}
    ) |>
    add_mean_bar() |>
    add_data_labels(label = label_vec, label_position = "above") |>
    adjust_colors(new_color_scheme(my_color)) |>
    sort_x_axis_labels(.reverse = TRUE) |>
    adjust_x_axis(rotate_labels = 90)
}

df_dxy <- read_csv(paths$dxy_csv)
df_dxy_long <- df_dxy |>
  pivot_longer(cols = -c(`...1`), names_to = "population_2", values_to = "dxy") |>
  rename(population1 = `...1`)

p_dxy <- format_heatmap(
  df_dxy_long,
  dxy,
  scales::label_number(accuracy = 1e-5, trim = TRUE)
)

df_fst_pair <- read_csv(paths$fst_pair)
df_fst_pair_long <- df_fst_pair |>
  pivot_longer(cols = -c(`...1`), names_to = "population_2", values_to = "fst") |>
  rename(population1 = `...1`)

p_fst_pair <- format_heatmap(df_fst_pair_long, fst)

df_pi <- read_csv(paths$pi_csv)
p_pi <- df_pi |>
  tidyplot(
    x = population,
    y = pi,
    color = pi
  ) |>
  add_sum_bar() |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_colors(
    new_colors = new_color_scheme(my_color),
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 1e-5, trim = TRUE)
  ) |>
  add_data_labels(label = round(df_pi$pi, 4), label_position = "above") |>
  sort_x_axis_labels(.reverse = TRUE)

df_tajima_d <- read_csv(paths$tajima_d_csv)
p_tajima_d <- format_bar(df_tajima_d, tajima_d, round(df_tajima_d$tajima_d, 4))

df_theta_w <- read_csv(paths$theta_w_csv)
p_df_theta_w <- format_bar(df_theta_w, theta_w, round(df_theta_w$theta_w, 4))

df_population_count <- read_csv(paths$population_count_csv)
p_df_population_count <- df_population_count |>
  tidyplot(
    x = group,
    y = count,
    color = count
  ) |>
  add_sum_bar() |>
  adjust_x_axis(rotate_labels = 90) |>
  adjust_colors(new_color_scheme(my_color)) |>
  add_data_labels(label = df_population_count$count, label_position = "above") |>
  sort_x_axis_labels(.reverse = TRUE)

p_all <- patchwork::wrap_plots(
  p_dxy,
  p_fst_pair,
  p_pi,
  p_tajima_d,
  p_df_theta_w,
  p_df_population_count
) +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  filename = out_file,
  plot = p_all,
  width = 30,
  height = 25,
  units = "cm"
)
