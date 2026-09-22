# Simulation S8: pointwise oracle-CSPD interval coverage.

s8_oracle_target <- function(grid, rho, alpha) {
  cutoff <- stats::qnorm(1 - alpha / 2)
  half_width <- sqrt(1 - rho^2) * cutoff - abs(rho) * max(abs(grid))
  if (half_width <= 0) stop("the selected grid has zero oracle common support.", call. = FALSE)
  gamma <- 2 * stats::pnorm(half_width) - 1
  second_moment <- 1 - 2 * half_width * stats::dnorm(half_width) / gamma
  data.frame(z = grid, truth = sin(grid) + 0.5 * second_moment, gamma = gamma)
}

run_s8 <- function(
  n = 400L,
  repetitions = 500L,
  rho = 0.7,
  alpha = 0.05,
  grid = seq(-0.75, 0.75, length.out = 31L),
  seed = 20260916L
) {
  grid <- assert_grid(grid)
  truth <- s8_oracle_target(grid, rho, alpha)
  is_valid <- s1_validity_rule(rho, alpha)
  set.seed(seed)
  intervals <- lapply(seq_len(repetitions), function(replication) {
    reference_data <- simulate_s1_data(n, rho)
    estimate_cspd_inference(reference_data, "x1", grid, is_valid, s1_prediction)
  })
  covered <- do.call(cbind, lapply(intervals, function(interval) {
    interval$lower <= truth$truth & truth$truth <= interval$upper
  }))
  estimates <- do.call(cbind, lapply(intervals, `[[`, "cspd"))
  standard_errors <- do.call(cbind, lapply(intervals, `[[`, "std_error"))
  gamma <- vapply(intervals, function(interval) unique(interval$gamma), numeric(1))
  data.frame(
    z = grid,
    truth = truth$truth,
    coverage = rowMeans(covered, na.rm = TRUE),
    mean_estimate = rowMeans(estimates, na.rm = TRUE),
    mean_std_error = rowMeans(standard_errors, na.rm = TRUE),
    oracle_gamma = truth$gamma,
    mean_sample_gamma = mean(gamma),
    repetitions = repetitions
  )
}

plot_s8 <- function(simulation, path = file.path("results", "s8-oracle-inference.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  graphics::plot(simulation$z, simulation$mean_estimate, type = "l", lwd = 2,
                 col = "#0072B2", xlab = expression(z), ylab = "CSPD")
  graphics::lines(simulation$z, simulation$truth, lwd = 2, lty = 2, col = "grey20")
  graphics::legend("bottomright", legend = c("mean estimate", "oracle target"),
                   col = c("#0072B2", "grey20"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  graphics::plot(simulation$z, simulation$coverage, type = "b", pch = 19, ylim = c(0.8, 1),
                 lwd = 2, col = "#009E73", xlab = expression(z), ylab = "Empirical coverage")
  graphics::abline(h = 0.95, lty = 2, lwd = 2, col = "grey20")
  path
}
