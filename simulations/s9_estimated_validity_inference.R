# Simulation S9: sample-split inference with an estimated validity domain.

s9_truncated_second_moment <- function(lower, upper) {
  probability <- stats::pnorm(upper) - stats::pnorm(lower)
  if (!is.finite(probability) || probability <= 0) return(NA_real_)
  x_phi <- function(x) ifelse(is.finite(x), x * stats::dnorm(x), 0)
  1 + (x_phi(lower) - x_phi(upper)) / probability
}

s9_target_from_interval <- function(grid, lower, upper) {
  second_moment <- s9_truncated_second_moment(lower, upper)
  data.frame(z = grid, target = sin(grid) + 0.5 * second_moment,
             gamma = stats::pnorm(upper) - stats::pnorm(lower))
}

s9_plugin_target <- function(conformal, grid) {
  coefficients <- stats::coef(conformal$model)
  intercept <- unname(coefficients[[1]])
  slope <- unname(coefficients[[2]])
  half_width <- conformal$radius - max(abs(grid))
  if (!is.finite(slope) || abs(slope) < 1e-8 || half_width <= 0) {
    return(s9_target_from_interval(grid, 1, 0))
  }
  bounds <- sort(c((-half_width - intercept) / slope, (half_width - intercept) / slope))
  s9_target_from_interval(grid, bounds[[1]], bounds[[2]])
}

s9_oracle_target <- function(grid, rho, alpha) {
  half_width <- sqrt(1 - rho^2) * stats::qnorm(1 - alpha / 2) - max(abs(grid))
  if (half_width <= 0) stop("the selected grid has zero oracle common support.", call. = FALSE)
  bounds <- c(-half_width / abs(rho), half_width / abs(rho))
  s9_target_from_interval(grid, bounds[[1]], bounds[[2]])
}

run_s9 <- function(
  n_training = 100L,
  n_calibration = 100L,
  n_evaluation = 400L,
  repetitions = 400L,
  rho = 0.7,
  alpha = 0.05,
  grid = seq(-0.75, 0.75, length.out = 31L),
  seed = 20260917L
) {
  grid <- assert_grid(grid)
  oracle_target <- s9_oracle_target(grid, rho, alpha)
  set.seed(seed)
  runs <- lapply(seq_len(repetitions), function(replication) {
    training_data <- simulate_s1_data(n_training, rho)
    calibration_data <- simulate_s1_data(n_calibration, rho)
    evaluation_data <- simulate_s1_data(n_evaluation, rho)
    conformal <- fit_split_conformal_validity(training_data, calibration_data, "x1", alpha)
    estimated_validity <- function(newdata) predict_validity(conformal, newdata)
    plugin_target <- s9_plugin_target(conformal, grid)
    interval <- estimate_cspd_inference(evaluation_data, "x1", grid, estimated_validity, s1_prediction)
    data.frame(
      z = grid,
      covered_plugin = interval$lower <= plugin_target$target & plugin_target$target <= interval$upper,
      covered_oracle = interval$lower <= oracle_target$target & oracle_target$target <= interval$upper,
      target_gap = plugin_target$target - oracle_target$target,
      gamma_plugin = plugin_target$gamma,
      gamma_evaluation = interval$gamma,
      radius = conformal$radius
    )
  })
  raw <- do.call(rbind, runs)
  aggregate(raw[, c("covered_plugin", "covered_oracle", "target_gap", "gamma_plugin", "gamma_evaluation")],
            by = list(z = raw$z), FUN = mean)
}

plot_s9 <- function(simulation, path = file.path("results", "s9-estimated-validity-inference.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  graphics::plot(simulation$z, simulation$covered_plugin, type = "b", pch = 19,
                 ylim = c(0, 1), lwd = 2, col = "#0072B2", xlab = expression(z),
                 ylab = "Empirical coverage")
  graphics::lines(simulation$z, simulation$covered_oracle, type = "b", pch = 17,
                 lwd = 2, lty = 2, col = "#D55E00")
  graphics::abline(h = 0.95, lwd = 2, lty = 3, col = "grey20")
  graphics::legend("bottom", legend = c("plug-in target", "oracle target", "nominal 95%"),
                   col = c("#0072B2", "#D55E00", "grey20"), lty = c(1, 2, 3),
                   pch = c(19, 17, NA), lwd = 2, bty = "n", cex = 0.8)
  graphics::plot(simulation$z, simulation$target_gap, type = "l", lwd = 2, col = "#D55E00",
                 xlab = expression(z), ylab = "Plug-in target minus oracle target")
  graphics::abline(h = 0, lwd = 2, lty = 2, col = "grey20")
  path
}
