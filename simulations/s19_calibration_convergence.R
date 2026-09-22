# Simulation S19: finite-sample convergence of learned residual validity to oracle validity.

s19_one <- function(split_size, repetition, n_evaluation, rho, alpha,
                    grid_half_width, grid_points, seed) {
  set.seed(seed + 100L * as.integer(split_size) + repetition)
  training <- simulate_s1_data(split_size, rho)
  calibration <- simulate_s1_data(split_size, rho)
  evaluation <- simulate_s1_data(n_evaluation, rho)
  grid <- seq(-grid_half_width, grid_half_width, length.out = grid_points)
  oracle <- s14_oracle_validity(rho, alpha)
  learned_fit <- fit_split_conformal_validity(training, calibration, "x1", alpha)
  learned <- function(newdata) predict_validity(learned_fit, newdata)
  oracle_h <- validity_matrix(evaluation, "x1", grid, oracle)
  learned_h <- validity_matrix(evaluation, "x1", grid, learned)
  data.frame(
    split_size = split_size, repetition = repetition,
    oracle_gamma = mean(rowSums(oracle_h) == ncol(oracle_h)),
    learned_gamma = mean(rowSums(learned_h) == ncol(learned_h)),
    oracle_mean_isc = mean(oracle_h), learned_mean_isc = mean(learned_h),
    radius_ratio = learned_fit$radius /
      (sqrt(1 - rho^2) * stats::qnorm(1 - alpha / 2))
  )
}

run_s19 <- function(split_sizes = c(100L, 200L, 500L, 1000L), repetitions = 200L,
                    n_evaluation = 2000L, rho = 0.8, alpha = 0.05,
                    grid_half_width = 1.1, grid_points = 21L,
                    seed = 20260924L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  design <- expand.grid(split_size = split_sizes, repetition = seq_len(repetitions))
  evaluate <- function(i) s19_one(design$split_size[[i]], design$repetition[[i]],
                                  n_evaluation, rho, alpha, grid_half_width,
                                  grid_points, seed)
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(design)), evaluate, mc.cores = workers)
  } else lapply(seq_len(nrow(design)), evaluate)
  do.call(rbind, pieces)
}

summarize_s19 <- function(simulation) {
  measures <- c("oracle_gamma", "learned_gamma", "oracle_mean_isc",
                "learned_mean_isc", "radius_ratio")
  pieces <- split(simulation, simulation$split_size)
  do.call(rbind, lapply(pieces, function(part) {
    output <- data.frame(split_size = part$split_size[[1]])
    for (measure in measures) {
      output[[measure]] <- mean(part[[measure]])
      output[[paste0(measure, "_mc_se")]] <- stats::sd(part[[measure]]) / sqrt(nrow(part))
    }
    output
  }))
}

plot_s19 <- function(simulation, path = file.path("results", "s19-calibration-convergence.png")) {
  summary <- summarize_s19(simulation)
  mass <- rbind(
    data.frame(split_size = summary$split_size, common_mass = summary$learned_gamma,
               rule = "Learned"),
    data.frame(split_size = summary$split_size, common_mass = summary$oracle_gamma,
               rule = "Oracle")
  )
  p_mass <- ggplot2::ggplot(
    mass, ggplot2::aes(split_size, common_mass, color = rule,
                       shape = rule, linetype = rule)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.2, stroke = 0.7) +
    ggplot2::scale_x_log10(breaks = summary$split_size) +
    ggplot2::scale_color_manual(values = c("Learned" = pdp_palette[["blue"]],
                                           "Oracle" = pdp_palette[["ink"]])) +
    ggplot2::scale_shape_manual(values = c("Learned" = 16, "Oracle" = 1)) +
    ggplot2::scale_linetype_manual(values = c("Learned" = "solid", "Oracle" = "dashed")) +
    ggplot2::labs(x = "Training and calibration size", y = "Common mass",
                  color = NULL, shape = NULL, linetype = NULL) +
    theme_pdp()
  p_radius <- ggplot2::ggplot(summary, ggplot2::aes(split_size, radius_ratio)) +
    ggplot2::geom_hline(yintercept = 1, color = pdp_palette[["ink"]], linetype = "dashed") +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.7) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.2) +
    ggplot2::scale_x_log10(breaks = summary$split_size) +
    ggplot2::labs(x = "Training and calibration size", y = "Learned/oracle radius") +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  save_pdp_panel_set(list(p_mass, p_radius), path,
                     widths = 5.2, heights = 4.0)
}
