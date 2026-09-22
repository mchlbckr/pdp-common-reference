# Simulation S5: compare oracle and split-conformal validity rules.

s5_oracle_validity <- function(rho, alpha = 0.05) {
  cutoff <- stats::qnorm(1 - alpha / 2)
  function(newdata) {
    residual <- (newdata$x1 - rho * newdata$x2) / sqrt(1 - rho^2)
    abs(residual) <= cutoff
  }
}

run_s5 <- function(n = 6000L, rho = 0.7, grid = seq(-0.75, 0.75, length.out = 41L), alpha = 0.05, seed = 20260912L) {
  set.seed(seed)
  data <- simulate_s1_data(n = n, rho = rho)
  split <- sample(rep(seq_len(3L), length.out = n))
  training_data <- data[split == 1L, , drop = FALSE]
  calibration_data <- data[split == 2L, , drop = FALSE]
  evaluation_data <- data[split == 3L, , drop = FALSE]

  conformal <- fit_split_conformal_validity(training_data, calibration_data, "x1", alpha)
  estimated_validity <- function(newdata) predict_validity(conformal, newdata)
  oracle_validity <- s5_oracle_validity(rho, alpha)

  isc_oracle <- estimate_isc(evaluation_data, "x1", grid, oracle_validity)
  isc_estimated <- estimate_isc(evaluation_data, "x1", grid, estimated_validity)
  cspd_oracle <- estimate_cspd(evaluation_data, "x1", grid, oracle_validity, s1_prediction)
  cspd_estimated <- estimate_cspd(evaluation_data, "x1", grid, estimated_validity, s1_prediction)

  data.frame(
    z = grid,
    isc_oracle = isc_oracle$isc,
    isc_estimated = isc_estimated$isc,
    cspd_oracle = cspd_oracle$cspd,
    cspd_estimated = cspd_estimated$cspd,
    gamma_oracle = cspd_oracle$gamma,
    gamma_estimated = cspd_estimated$gamma,
    radius = conformal$radius,
    rho = rho
  )
}
