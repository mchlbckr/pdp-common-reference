# Simulation S22: empirical RCPD convergence to a population approximation.

s22_oracle_target <- function(n_oracle = 100000L, rho = 0.8, alpha = 0.05,
                              epsilon = 0.10,
                              grid = seq(-1.1, 1.1, length.out = 21L),
                              seed = 20260925L) {
  population <- simulate_s1_data(n_oracle, rho, seed = seed)
  valid <- s14_oracle_validity(rho, alpha)
  estimate_relaxed_cspd(
    population, "x1", grid, valid, s18_prediction,
    epsilon = epsilon, tolerance = 2e-5
  )$relaxed_cspd
}

run_s22 <- function(sample_sizes = c(250L, 500L, 1000L, 2000L, 4000L),
                    repetitions = 100L, rho = 0.8, alpha = 0.05,
                    epsilon = 0.10,
                    grid = seq(-1.1, 1.1, length.out = 21L),
                    n_oracle = 100000L, seed = 20260925L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  oracle <- s22_oracle_target(n_oracle, rho, alpha, epsilon, grid, seed)
  design <- expand.grid(
    sample_size = sample_sizes, repetition = seq_len(repetitions),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  valid <- s14_oracle_validity(rho, alpha)
  evaluate <- function(index) {
    setting <- design[index, ]
    data <- simulate_s1_data(setting$sample_size, rho, seed = seed + index)
    estimate <- estimate_relaxed_cspd(
      data, "x1", grid, valid, s18_prediction,
      epsilon = epsilon, tolerance = 2e-5
    )
    data.frame(
      sample_size = setting$sample_size, repetition = setting$repetition,
      sup_norm_error = max(abs(estimate$relaxed_cspd - oracle)),
      active_constraints = unique(estimate$active_constraints),
      binding_constraints = unique(estimate$binding_constraints),
      strict_complementarity = unique(estimate$strict_complementarity),
      stringsAsFactors = FALSE
    )
  }
  indices <- seq_len(nrow(design))
  output <- if (workers > 1L && .Platform$OS.type != "windows") {
    do.call(rbind, parallel::mclapply(indices, evaluate, mc.cores = workers))
  } else {
    do.call(rbind, lapply(indices, evaluate))
  }
  attr(output, "oracle_curve") <- oracle
  attr(output, "grid") <- grid
  output
}

summarize_s22 <- function(simulation) {
  pieces <- split(simulation, simulation$sample_size)
  summary <- do.call(rbind, lapply(pieces, function(part) {
    data.frame(
      sample_size = part$sample_size[[1]],
      mean_sup_norm_error = mean(part$sup_norm_error),
      error_mc_se = stats::sd(part$sup_norm_error) / sqrt(nrow(part)),
      strict_complementarity_rate = mean(part$strict_complementarity),
      mean_active_constraints = mean(part$active_constraints),
      stringsAsFactors = FALSE
    )
  }))
  slope <- unname(stats::coef(stats::lm(
    log(mean_sup_norm_error) ~ log(sample_size), data = summary
  ))[[2]])
  summary$log_log_slope <- slope
  summary
}

plot_s22 <- function(summary, path = file.path("results", "s22-rcpd-rate.png")) {
  reference <- data.frame(
    sample_size = range(summary$sample_size),
    error = summary$mean_sup_norm_error[[1]] *
      (range(summary$sample_size) / min(summary$sample_size))^(-0.5)
  )
  figure <- ggplot2::ggplot(
    summary, ggplot2::aes(sample_size, mean_sup_norm_error)
  ) +
    ggplot2::geom_line(data = reference, ggplot2::aes(sample_size, error),
                       inherit.aes = FALSE, linetype = "dashed",
                       color = pdp_palette[["mid_grey"]], linewidth = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_sup_norm_error - 2 * error_mc_se,
                   ymax = mean_sup_norm_error + 2 * error_mc_se),
      width = 0.04, linewidth = 0.6, color = pdp_palette[["blue"]]
    ) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.75) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.2) +
    ggplot2::scale_x_log10(breaks = summary$sample_size) +
    ggplot2::scale_y_log10() +
    ggplot2::annotate(
      "text", x = max(summary$sample_size), y = max(summary$mean_sup_norm_error),
      label = sprintf("estimated slope %.2f", unique(summary$log_log_slope)),
      hjust = 1, vjust = 1.2, size = 3
    ) +
    ggplot2::labs(x = "Evaluation sample size (log scale)",
                  y = "Mean sup-norm error (log scale)") +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  save_pdp_plot(figure, path, width = 7.2, height = 4.8)
}
