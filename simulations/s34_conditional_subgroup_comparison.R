# Simulation S34: a transparent conditional-subgroup PDP comparison for the
# additive composition design. The comparator follows the cs-PDP principle of
# Molnar et al. (2024): partition the context space so that the focal feature is
# more homogeneous within groups, then report one PDP per group. It is not
# collapsed into a single global curve because that would change its estimand.

run_s34 <- function(n_evaluation = 1000L, sigma = 0.45, alpha = 0.05,
                    grid = seq(-0.75, 0.75, length.out = 21L), beta = 1.5,
                    strength = 4, seed = 20260926L, maxdepth = 2L,
                    minbucket_fraction = 0.12) {
  grid <- assert_grid(grid)
  data <- simulate_s23_data(n_evaluation, sigma, seed)
  tree <- rpart::rpart(
    x ~ u, data = data, method = "anova",
    control = rpart::rpart.control(
      cp = 0, maxdepth = maxdepth,
      minbucket = max(20L, floor(minbucket_fraction * n_evaluation)),
      minsplit = max(40L, floor(2 * minbucket_fraction * n_evaluation)),
      xval = 0
    )
  )
  raw_group <- tree$where
  group_order <- names(sort(tapply(data$u, raw_group, mean)))
  group_index <- match(as.character(raw_group), group_order)
  data$group <- factor(
    sprintf("G%d", group_index),
    levels = sprintf("G%d", seq_along(group_order))
  )
  prediction <- s23_prediction(
    direction = 1, sigma = sigma, alpha = alpha,
    beta = beta, strength = strength
  )

  curve_parts <- lapply(levels(data$group), function(group_label) {
    selected <- data$group == group_label
    weights <- as.numeric(selected) / sum(selected)
    curve <- estimate_weighted_pdp(
      data, "x", grid, prediction, weights
    )$weighted_pdp
    focal_range <- stats::quantile(
      data$x[selected], probs = c(0.25, 0.75), names = FALSE, type = 8
    )
    data.frame(
      group = group_label,
      z = grid,
      value = curve - mean(curve),
      dense = grid >= focal_range[[1]] & grid <= focal_range[[2]],
      stringsAsFactors = FALSE
    )
  })
  curves <- do.call(rbind, curve_parts)
  group_summary <- do.call(rbind, lapply(levels(data$group), function(group_label) {
    selected <- data$group == group_label
    focal_quantiles <- stats::quantile(
      data$x[selected], probs = c(0.25, 0.50, 0.75),
      names = FALSE, type = 8
    )
    data.frame(
      group = group_label,
      n = sum(selected),
      fraction = mean(selected),
      mean_u = mean(data$u[selected]),
      x_q25 = focal_quantiles[[1]],
      x_median = focal_quantiles[[2]],
      x_q75 = focal_quantiles[[3]],
      stringsAsFactors = FALSE
    )
  }))
  list(
    data = data,
    curves = curves,
    groups = group_summary,
    truth = data.frame(z = grid, value = s23_center(s23_truth(grid))),
    settings = list(
      n_evaluation = n_evaluation, sigma = sigma, alpha = alpha,
      grid = grid, beta = beta, strength = strength, seed = seed,
      maxdepth = maxdepth, minbucket_fraction = minbucket_fraction
    )
  )
}

summarize_s34 <- function(simulation) simulation$groups

plot_s34 <- function(simulation, directory = "results") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  group_levels <- levels(simulation$data$group)
  colors <- setNames(
    c(
      pdp_palette[["blue"]], pdp_palette[["orange"]],
      pdp_palette[["green"]], pdp_palette[["purple"]]
    )[seq_along(group_levels)],
    group_levels
  )

  display_data <- simulation$data
  if (nrow(display_data) > 600L) {
    set.seed(simulation$settings$seed + 1L)
    display_data <- display_data[sample.int(nrow(display_data), 600L), ]
  }
  panel_partition <- ggplot2::ggplot(
    display_data, ggplot2::aes(u, x, color = group)
  ) +
    ggplot2::geom_point(size = 1.0, alpha = 0.48) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = "Context U", y = "Focal feature X", color = "CART subgroup"
    ) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  full_curves <- simulation$curves
  dense_curves <- simulation$curves[simulation$curves$dense, , drop = FALSE]
  panel_curves <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = full_curves,
      ggplot2::aes(z, value, color = group, group = group),
      linewidth = 0.55, alpha = 0.30
    ) +
    ggplot2::geom_line(
      data = dense_curves,
      ggplot2::aes(z, value, color = group, group = group),
      linewidth = 1.0
    ) +
    ggplot2::geom_line(
      data = simulation$truth,
      ggplot2::aes(z, value),
      color = pdp_palette[["ink"]], linetype = "dashed", linewidth = 0.75
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = expression(z), y = "Centered subgroup PD", color = "CART subgroup"
    ) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  paths <- file.path(
    directory,
    c(
      "s34-conditional-subgroup-comparison-a.png",
      "s34-conditional-subgroup-comparison-b.png"
    )
  )
  save_pdp_plot(panel_partition, paths[[1]], width = 5.2, height = 3.8)
  save_pdp_plot(panel_curves, paths[[2]], width = 5.2, height = 3.8)
  paths
}
