# Simulation S25: separate evaluation variation from validity-learning error.

s25_one <- function(split_size, repetition, rho, alpha, grid, seed) {
  set.seed(seed + 100000L * as.integer(split_size) + repetition)
  training <- simulate_s1_data(split_size, rho)
  calibration <- simulate_s1_data(split_size, rho)
  evaluation <- simulate_s1_data(4L * split_size, rho)
  fit <- fit_split_conformal_validity(training, calibration, "x1", alpha)
  learned_validity <- function(newdata) predict_validity(fit, newdata)
  plugin_target <- s9_plugin_target(fit, grid)
  oracle_target <- s9_oracle_target(grid, rho, alpha)
  estimate <- estimate_cspd_inference(
    evaluation, "x1", grid, learned_validity, s1_prediction
  )
  data.frame(
    split_size = split_size,
    total_n = 6L * split_size,
    repetition = repetition,
    z = grid,
    estimate = estimate$cspd,
    plugin_target = plugin_target$target,
    oracle_target = oracle_target$target,
    standard_error = estimate$std_error,
    evaluation_error = estimate$cspd - plugin_target$target,
    validity_error = plugin_target$target - oracle_target$target,
    total_error = estimate$cspd - oracle_target$target,
    covered_plugin = estimate$lower <= plugin_target$target &
      plugin_target$target <= estimate$upper,
    covered_oracle = estimate$lower <= oracle_target$target &
      oracle_target$target <= estimate$upper,
    plugin_gamma = plugin_target$gamma,
    evaluation_gamma = estimate$gamma,
    radius_ratio = fit$radius /
      (sqrt(1 - rho^2) * stats::qnorm(1 - alpha / 2)),
    stringsAsFactors = FALSE
  )
}

run_s25 <- function(split_sizes = c(75L, 150L, 300L, 600L, 1200L),
                    repetitions = 200L, rho = 0.7, alpha = 0.05,
                    grid = seq(-0.75, 0.75, length.out = 11L),
                    seed = 20260927L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  grid <- assert_grid(grid)
  design <- expand.grid(split_size = split_sizes,
                        repetition = seq_len(repetitions))
  evaluate <- function(index) s25_one(
    design$split_size[[index]], design$repetition[[index]], rho, alpha,
    grid, seed
  )
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(design)), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(nrow(design)), evaluate)
  }
  do.call(rbind, pieces)
}

summarize_s25 <- function(simulation) {
  pieces <- split(simulation, simulation$split_size)
  do.call(rbind, lapply(pieces, function(part) {
    valid <- is.finite(part$evaluation_error) & is.finite(part$validity_error) &
      is.finite(part$total_error)
    part <- part[valid, , drop = FALSE]
    data.frame(
      split_size = part$split_size[[1]], total_n = part$total_n[[1]],
      evaluation_rmse = sqrt(mean(part$evaluation_error^2)),
      validity_rmse = sqrt(mean(part$validity_error^2)),
      total_rmse = sqrt(mean(part$total_error^2)),
      mean_validity_bias = mean(part$validity_error),
      mean_abs_validity_bias = mean(abs(part$validity_error)),
      empirical_evaluation_sd = stats::sd(part$evaluation_error),
      mean_reported_se = mean(part$standard_error),
      plugin_coverage = mean(part$covered_plugin),
      oracle_coverage = mean(part$covered_oracle),
      radius_ratio = mean(part$radius_ratio),
      plugin_gamma = mean(part$plugin_gamma),
      evaluation_gamma = mean(part$evaluation_gamma),
      usable_fraction = nrow(part) / nrow(simulation[simulation$split_size ==
                                                       part$split_size[[1]], ]),
      stringsAsFactors = FALSE
    )
  }))
}

plot_s25 <- function(summary,
                     path = file.path("results", "s25-error-decomposition.png")) {
  errors <- rbind(
    data.frame(total_n = summary$total_n, value = summary$evaluation_rmse,
               component = "Evaluation RMSE"),
    data.frame(total_n = summary$total_n, value = summary$validity_rmse,
               component = "Validity-target RMSE"),
    data.frame(total_n = summary$total_n, value = summary$total_rmse,
               component = "Total oracle-target RMSE")
  )
  errors$component <- factor(
    errors$component,
    levels = c("Evaluation RMSE", "Validity-target RMSE",
               "Total oracle-target RMSE")
  )
  colors <- c(
    `Evaluation RMSE` = pdp_palette[["blue"]],
    `Validity-target RMSE` = pdp_palette[["orange"]],
    `Total oracle-target RMSE` = pdp_palette[["ink"]]
  )
  p_error <- ggplot2::ggplot(
    errors, ggplot2::aes(total_n, value, color = component,
                         shape = component, linetype = component)
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.1) +
    ggplot2::scale_x_log10(breaks = summary$total_n) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_shape_manual(values = c(16, 17, 1)) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
    ggplot2::labs(x = "Total split-sample size", y = "Root mean squared error",
                  color = NULL, shape = NULL, linetype = NULL) +
    theme_pdp()

  variance <- rbind(
    data.frame(total_n = summary$total_n,
               value = summary$empirical_evaluation_sd,
               quantity = "Empirical evaluation SD"),
    data.frame(total_n = summary$total_n,
               value = summary$mean_reported_se,
               quantity = "Mean reported SE")
  )
  p_variance <- ggplot2::ggplot(
    variance, ggplot2::aes(total_n, value, color = quantity,
                           shape = quantity, linetype = quantity)
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.1) +
    ggplot2::scale_x_log10(breaks = summary$total_n) +
    ggplot2::scale_color_manual(values = c(
      `Empirical evaluation SD` = pdp_palette[["ink"]],
      `Mean reported SE` = pdp_palette[["sky"]]
    )) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
    ggplot2::labs(x = "Total split-sample size",
                  y = "Evaluation standard deviation",
                  color = NULL, shape = NULL, linetype = NULL) +
    theme_pdp()

  coverage <- rbind(
    data.frame(total_n = summary$total_n, coverage = summary$plugin_coverage,
               target = "Learned plug-in target"),
    data.frame(total_n = summary$total_n, coverage = summary$oracle_coverage,
               target = "Fixed oracle target")
  )
  p_coverage <- ggplot2::ggplot(
    coverage, ggplot2::aes(total_n, coverage, color = target,
                           shape = target, linetype = target)
  ) +
    ggplot2::geom_hline(yintercept = 0.95, color = pdp_palette[["ink"]],
                        linetype = "dotted", linewidth = 0.55) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.1) +
    ggplot2::scale_x_log10(breaks = summary$total_n) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_color_manual(values = c(
      `Learned plug-in target` = pdp_palette[["blue"]],
      `Fixed oracle target` = pdp_palette[["orange"]]
    )) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
    ggplot2::labs(x = "Total split-sample size", y = "Pointwise coverage",
                  color = NULL, shape = NULL, linetype = NULL) +
    theme_pdp()

  save_pdp_panel_set(
    list(p_error, p_variance, p_coverage), path,
    widths = 5.2, heights = 4.0
  )
}
