# Simulation S10: outer bootstrap for the fixed oracle CSPD target.

s10_estimate_cspd <- function(training_data, calibration_data, evaluation_data, grid, alpha) {
  conformal <- fit_split_conformal_validity(training_data, calibration_data, "x1", alpha)
  estimated_validity <- function(newdata) predict_validity(conformal, newdata)
  estimate_cspd_inference(evaluation_data, "x1", grid, estimated_validity, s1_prediction)
}

s10_resample_rows <- function(data) {
  data[sample.int(nrow(data), nrow(data), replace = TRUE), , drop = FALSE]
}

run_s10 <- function(
  n_training = 100L,
  n_calibration = 100L,
  n_evaluation = 400L,
  repetitions = 100L,
  bootstrap_repetitions = 60L,
  rho = 0.7,
  alpha = 0.05,
  grid = seq(-0.75, 0.75, length.out = 21L),
  seed = 20260918L
) {
  grid <- assert_grid(grid)
  oracle_target <- s9_oracle_target(grid, rho, alpha)$target
  set.seed(seed)
  runs <- lapply(seq_len(repetitions), function(replication) {
    training_data <- simulate_s1_data(n_training, rho)
    calibration_data <- simulate_s1_data(n_calibration, rho)
    evaluation_data <- simulate_s1_data(n_evaluation, rho)
    point <- s10_estimate_cspd(training_data, calibration_data, evaluation_data, grid, alpha)
    bootstrap <- replicate(bootstrap_repetitions, {
      bootstrap_estimate <- s10_estimate_cspd(
        s10_resample_rows(training_data), s10_resample_rows(calibration_data),
        s10_resample_rows(evaluation_data), grid, alpha
      )
      bootstrap_estimate$cspd
    })
    lower <- apply(bootstrap, 1, stats::quantile, probs = 0.025, na.rm = TRUE, names = FALSE)
    upper <- apply(bootstrap, 1, stats::quantile, probs = 0.975, na.rm = TRUE, names = FALSE)
    data.frame(
      z = grid,
      covered_wald = point$lower <= oracle_target & oracle_target <= point$upper,
      covered_bootstrap = lower <= oracle_target & oracle_target <= upper,
      bootstrap_width = upper - lower
    )
  })
  raw <- do.call(rbind, runs)
  aggregate(raw[, c("covered_wald", "covered_bootstrap", "bootstrap_width")],
            by = list(z = raw$z), FUN = mean)
}

plot_s10 <- function(simulation, path = file.path("results", "s10-outer-bootstrap.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  graphics::plot(simulation$z, simulation$covered_wald, type = "b", pch = 17,
                 ylim = c(0, 1), lwd = 2, lty = 2, col = "#D55E00", xlab = expression(z),
                 ylab = "Oracle-target coverage")
  graphics::lines(simulation$z, simulation$covered_bootstrap, type = "b", pch = 19,
                 lwd = 2, col = "#0072B2")
  graphics::abline(h = 0.95, lwd = 2, lty = 3, col = "grey20")
  graphics::legend("bottom", legend = c("plug-in Wald", "outer bootstrap", "nominal 95%"),
                   col = c("#D55E00", "#0072B2", "grey20"), lty = c(2, 1, 3),
                   pch = c(17, 19, NA), lwd = 2, bty = "n", cex = 0.8)
  graphics::plot(simulation$z, simulation$bootstrap_width, type = "l", lwd = 2, col = "#0072B2",
                 xlab = expression(z), ylab = "Mean percentile-interval width")
  path
}
