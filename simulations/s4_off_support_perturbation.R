# Simulation S4: models agree on the validity domain but differ outside it.

simulate_s4_data <- function(n, rho, seed = NULL) {
  simulate_s1_data(n = n, rho = rho, seed = seed)
}

s4_base_model <- function(newdata) {
  stats::plogis(0.8 * newdata$x1 - 0.4 * newdata$x2 + 0.25 * newdata$x1 * newdata$x2)
}

s4_perturbed_model <- function(newdata, rho, alpha = 0.05, strength = 4) {
  cutoff <- stats::qnorm(1 - alpha / 2)
  residual <- (newdata$x2 - rho * newdata$x1) / sqrt(1 - rho^2)
  excess <- pmax(abs(residual) - cutoff, 0)
  stats::plogis(stats::qlogis(s4_base_model(newdata)) + strength * excess)
}

run_s4 <- function(
  n = 5000L,
  rho = 0.7,
  grid = seq(-0.75, 0.75, length.out = 41L),
  alpha = 0.05,
  strength = 4,
  seed = 20260911L
) {
  grid <- assert_grid(grid)
  reference_data <- simulate_s4_data(n = n, rho = rho, seed = seed)
  is_valid <- s1_validity_rule(rho = rho, alpha = alpha)
  base_prediction <- s4_base_model
  perturbed_prediction <- function(newdata) {
    s4_perturbed_model(newdata, rho = rho, alpha = alpha, strength = strength)
  }

  pdp_base <- estimate_pdp(reference_data, "x1", grid, base_prediction)
  pdp_perturbed <- estimate_pdp(reference_data, "x1", grid, perturbed_prediction)
  isc <- estimate_isc(reference_data, "x1", grid, is_valid)
  cspd_base <- estimate_cspd(reference_data, "x1", grid, is_valid, base_prediction)
  cspd_perturbed <- estimate_cspd(reference_data, "x1", grid, is_valid, perturbed_prediction)

  empirical_bound <- vapply(grid, function(z) {
    query_data <- make_queries(reference_data, "x1", z)
    invalid <- !is_valid(query_data)
    if (!any(invalid)) return(0)
    isc_row <- isc[isc$z == z, , drop = FALSE]
    isc_row$per * max(abs(base_prediction(query_data)[invalid] - perturbed_prediction(query_data)[invalid]))
  }, numeric(1))

  data.frame(
    z = grid,
    pdp_base = pdp_base$pdp,
    pdp_perturbed = pdp_perturbed$pdp,
    pdp_gap = abs(pdp_base$pdp - pdp_perturbed$pdp),
    bound = empirical_bound,
    isc = isc$isc,
    per = isc$per,
    cspd_base = cspd_base$cspd,
    cspd_perturbed = cspd_perturbed$cspd,
    cspd_gap = abs(cspd_base$cspd - cspd_perturbed$cspd),
    gamma = cspd_base$gamma,
    n_common = cspd_base$n_common,
    rho = rho,
    alpha = alpha,
    strength = strength
  )
}

plot_s4 <- function(simulation, path = file.path("results", "s4-off-support-perturbation.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1800, height = 650, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  graphics::plot(
    simulation$z, simulation$pdp_base,
    type = "l", lwd = 2, col = "grey20",
    xlab = expression(z), ylab = "Average prediction"
  )
  graphics::lines(simulation$z, simulation$pdp_perturbed, lwd = 2, col = "#D55E00", lty = 2)
  graphics::legend(
    "bottomright", legend = c(expression(f), expression(tilde(f))),
    col = c("grey20", "#D55E00"), lty = c(1, 2), lwd = 2, bty = "n"
  )

  graphics::plot(
    simulation$z, simulation$pdp_gap,
    type = "l", ylim = range(c(0, simulation$bound)), lwd = 2, col = "#D55E00",
    xlab = expression(z), ylab = "Difference"
  )
  graphics::lines(simulation$z, simulation$bound, lwd = 2, col = "#0072B2", lty = 2)
  graphics::legend(
    "topleft", legend = c("|PDP(f)-PDP(f~)|", "PER × maximum off-support gap"),
    col = c("#D55E00", "#0072B2"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.72
  )

  graphics::plot(
    simulation$z, simulation$cspd_base,
    type = "l", lwd = 2, col = "grey20",
    xlab = expression(z), ylab = "Common-support prediction"
  )
  graphics::lines(simulation$z, simulation$cspd_perturbed, lwd = 2, col = "#D55E00", lty = 2)
  graphics::mtext(
    sprintf("gamma(G) = %.3f; max CSPD gap = %.2e", unique(simulation$gamma), max(simulation$cspd_gap, na.rm = TRUE)),
    side = 3, line = -1.6, cex = 0.72, adj = 0.98
  )

  path
}
