# Simulation S6: conceptual comparison with ALE and a conditional curve.

run_s6 <- function(
  n = 5000L,
  rho = 0.7,
  alpha = 0.05,
  grid = seq(-0.75, 0.75, length.out = 41L),
  ale_bins = 20L,
  seed = 20260915L
) {
  grid <- assert_grid(grid)
  reference_data <- simulate_s1_data(n = n, rho = rho, seed = seed)
  is_valid <- s1_validity_rule(rho, alpha)
  pdp <- estimate_pdp(reference_data, "x1", grid, s1_prediction)
  cspd <- estimate_cspd(reference_data, "x1", grid, is_valid, s1_prediction)
  isc <- estimate_isc(reference_data, "x1", grid, is_valid)

  # ALE uses the full observed x1 range, so that each local interval remains an
  # observational rather than an interventional comparison.
  breaks <- unique(stats::quantile(reference_data$x1, probs = seq(0, 1, length.out = ale_bins + 1L)))
  ale <- estimate_ale(reference_data, "x1", breaks, s1_prediction)

  # This is the exact conditional estimand for the S1 data-generating process;
  # it is included as an estimand comparator, not presented as an estimator.
  conditional <- sin(grid) + 0.5 * (rho^2 * grid^2 + 1 - rho^2)

  list(
    curves = list(pdp = pdp, cspd = cspd, ale = ale,
                  conditional = data.frame(z = grid, conditional = conditional)),
    diagnostics = cbind(isc, gamma = cspd$gamma),
    summary = data.frame(
      rho = rho,
      gamma = unique(cspd$gamma),
      minimum_isc = min(isc$isc),
      maximum_per = max(isc$per),
      stringsAsFactors = FALSE
    )
  )
}

plot_s6 <- function(simulation, path = file.path("results", "s6-baseline-comparison.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  curves <- simulation$curves
  diagnostics <- simulation$diagnostics
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

  yrange <- range(curves$pdp$pdp, curves$cspd$cspd, curves$ale$ale, curves$conditional$conditional)
  graphics::plot(curves$pdp$z, curves$pdp$pdp, type = "l", lwd = 2,
                 col = "grey20", ylim = yrange, xlab = expression(z),
                 ylab = "Feature-effect estimand")
  graphics::lines(curves$cspd$z, curves$cspd$cspd, lwd = 2, lty = 2, col = "#0072B2")
  graphics::lines(curves$ale$z, curves$ale$ale, lwd = 2, lty = 3, col = "#009E73")
  graphics::lines(curves$conditional$z, curves$conditional$conditional, lwd = 2, lty = 4, col = "#D55E00")
  graphics::legend("bottomright", legend = c("PDP", "CSPD", "ALE", "conditional"),
                   col = c("grey20", "#0072B2", "#009E73", "#D55E00"),
                   lty = c(1, 2, 3, 4), lwd = 2, bty = "n", cex = 0.8)

  graphics::plot(diagnostics$z, diagnostics$isc, type = "l", ylim = c(0, 1),
                 lwd = 2, col = "#009E73", xlab = expression(z),
                 ylab = "Support coverage")
  graphics::lines(diagnostics$z, 1 - diagnostics$isc, lwd = 2, lty = 2, col = "#D55E00")
  graphics::abline(h = unique(diagnostics$gamma), lwd = 2, lty = 3, col = "#0072B2")
  graphics::legend("bottom", legend = c("ISC", "PER", expression(gamma(G))),
                   col = c("#009E73", "#D55E00", "#0072B2"), lty = c(1, 2, 3),
                   lwd = 2, bty = "n", cex = 0.8)
  path
}
