# Simulation S3: common-reference effects under model interactions.

s3_prediction <- function(newdata, beta) {
  sin(newdata$x1) + 0.5 * newdata$x2^2 + beta * newdata$x1 * newdata$x2
}

run_s3 <- function(
  n = 4000L,
  rho = 0.7,
  betas = c(0, 0.75, 1.5),
  grid = seq(-0.75, 0.75, length.out = 41L),
  alpha = 0.05,
  seed = 20260914L
) {
  reference_data <- simulate_s1_data(n = n, rho = rho, seed = seed)
  is_valid <- s1_validity_rule(rho, alpha)
  common <- common_support_indicator(reference_data, "x1", grid, is_valid)
  common_mean_x2 <- mean(reference_data$x2[common])

  output <- lapply(betas, function(beta) {
    prediction <- function(newdata) s3_prediction(newdata, beta)
    pdp <- estimate_pdp(reference_data, "x1", grid, prediction)
    cspd <- estimate_cspd(reference_data, "x1", grid, is_valid, prediction)
    data.frame(
      beta = beta,
      z = grid,
      pdp = pdp$pdp,
      cspd = cspd$cspd,
      gamma = cspd$gamma,
      n_common = cspd$n_common,
      common_reference_mean_x2 = common_mean_x2,
      conditional_context_mean_x2 = rho * grid,
      conditional_effect = sin(grid) + 0.5 * ((rho * grid)^2 + 1 - rho^2) + beta * rho * grid^2
    )
  })
  do.call(rbind, output)
}

plot_s3 <- function(simulation, path = file.path("results", "s3-interactions.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  betas <- sort(unique(simulation$beta))
  grDevices::png(path, width = 1800, height = 1200, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (beta in betas) {
    s <- simulation[simulation$beta == beta, ]
    graphics::plot(s$z, s$pdp, type = "l", lwd = 2, col = "grey20", xlab = expression(z), ylab = "Feature effect")
    graphics::lines(s$z, s$cspd, lwd = 2, lty = 2, col = "#0072B2")
    graphics::lines(s$z, s$conditional_effect, lwd = 2, lty = 3, col = "#D55E00")
    graphics::legend("bottomright", legend = c("PDP", "CSPD", "conditional context"), col = c("grey20", "#0072B2", "#D55E00"), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.75)
  }
  s <- simulation[simulation$beta == betas[[1]], ]
  graphics::plot(s$z, s$conditional_context_mean_x2, type = "l", lwd = 2, col = "#D55E00", xlab = expression(z), ylab = expression("mean context " * X[2]))
  graphics::abline(h = unique(s$common_reference_mean_x2), col = "#0072B2", lwd = 2, lty = 2)
  graphics::legend("topleft", legend = c("conditional context mean", "CSPD common-reference mean"), col = c("#D55E00", "#0072B2"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.75)
  path
}
