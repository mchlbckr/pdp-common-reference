simulate_s2_data <- function(n, sigma = 0.12, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::runif(n, -2, 2)
  data.frame(x1 = x1, x2 = x1^2 + stats::rnorm(n, sd = sigma))
}

s2_prediction <- function(newdata) sin(newdata$x1) + 0.4 * newdata$x2

run_s2 <- function(n = 6000L, sigma = 0.12, grid = seq(-2, 2, length.out = 61L), alpha = 0.05, seed = 20260913L) {
  set.seed(seed)
  data <- simulate_s2_data(n, sigma)
  split <- sample(rep(1:3, length.out = n))
  train <- data[split == 1, , drop = FALSE]
  calibration <- data[split == 2, , drop = FALSE]
  evaluation <- data[split == 3, , drop = FALSE]
  oracle <- function(newdata) abs(newdata$x2 - newdata$x1^2) <= stats::qnorm(1 - alpha / 2) * sigma
  linear <- fit_split_conformal_validity(train, calibration, "x1", alpha)
  knn <- fit_knn_validity(train, alpha = alpha)
  linear_valid <- function(newdata) predict_validity(linear, newdata)
  knn_valid <- function(newdata) predict_knn_validity(knn, newdata)
  oracle_isc <- estimate_isc(evaluation, "x1", grid, oracle)
  linear_isc <- estimate_isc(evaluation, "x1", grid, linear_valid)
  knn_isc <- estimate_isc(evaluation, "x1", grid, knn_valid)
  data.frame(z = grid, isc_oracle = oracle_isc$isc, isc_linear = linear_isc$isc, isc_knn = knn_isc$isc)
}

plot_s2 <- function(simulation, path = file.path("results", "s2-nonlinear-manifold.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1400, height = 850, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(simulation$z, simulation$isc_oracle, type = "l", ylim = c(0, 1), lwd = 3,
    col = "grey20", xlab = expression(z), ylab = "Interventional Support Coverage")
  graphics::lines(simulation$z, simulation$isc_linear, lwd = 2, col = "#D55E00", lty = 2)
  graphics::lines(simulation$z, simulation$isc_knn, lwd = 2, col = "#0072B2", lty = 3)
  graphics::legend("bottom", legend = c("Oracle", "linear split conformal", "local kNN"),
    col = c("grey20", "#D55E00", "#0072B2"), lty = c(1, 2, 3), lwd = c(3, 2, 2), bty = "n")
  path
}
