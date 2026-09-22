# Simulation S1: additive model under controlled feature correlation.

simulate_s1_data <- function(n, rho, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("rho must lie strictly between -1 and 1.", call. = FALSE)
  }

  x1 <- stats::rnorm(n)
  x2 <- rho * x1 + sqrt(1 - rho^2) * stats::rnorm(n)
  data.frame(x1 = x1, x2 = x2)
}

s1_prediction <- function(newdata) {
  sin(newdata$x1) + 0.5 * newdata$x2^2
}

s1_validity_rule <- function(rho, alpha = 0.05) {
  cutoff <- stats::qnorm(1 - alpha / 2)
  function(newdata) {
    residual <- (newdata$x2 - rho * newdata$x1) / sqrt(1 - rho^2)
    abs(residual) <= cutoff
  }
}

run_s1 <- function(n = 2000L, rhos = c(0, 0.5, 0.9, 0.99), grid = seq(-2, 2, length.out = 41L), alpha = 0.05, seed = 20260910L) {
  grid <- assert_grid(grid)
  output <- lapply(seq_along(rhos), function(index) {
    rho <- rhos[[index]]
    reference_data <- simulate_s1_data(n = n, rho = rho, seed = seed + index)
    is_valid <- s1_validity_rule(rho = rho, alpha = alpha)

    pdp <- estimate_pdp(reference_data, "x1", grid, s1_prediction)
    isc <- estimate_isc(reference_data, "x1", grid, is_valid)
    cspd <- estimate_cspd(reference_data, "x1", grid, is_valid, s1_prediction)

    data.frame(
      rho = rho,
      z = grid,
      pdp = pdp$pdp,
      isc = isc$isc,
      per = isc$per,
      cspd = cspd$cspd,
      gamma = cspd$gamma,
      n_common = cspd$n_common
    )
  })

  do.call(rbind, output)
}
plot_s1 <- function(simulation, path = file.path("results", "s1-additive-correlation.png")) {
  required <- c("rho", "z", "pdp", "isc", "per", "cspd", "gamma")
  if (!all(required %in% names(simulation))) {
    stop("simulation does not have the required S1 columns.", call. = FALSE)
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  rho_values <- sort(unique(simulation$rho))
  grDevices::png(path, width = 1800, height = 1200, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (rho in rho_values) {
    subset <- simulation[simulation$rho == rho, ]
    graphics::plot(
      subset$z, subset$pdp,
      type = "l", lwd = 2, col = "grey20",
      xlab = expression(z), ylab = "Feature effect"
    )
    graphics::lines(subset$z, subset$cspd, lwd = 2, col = "#0072B2", lty = 2)
    graphics::abline(h = 0, col = "grey85")
    graphics::mtext(
      sprintf("Common-support mass: %.3f", unique(subset$gamma)),
      side = 3, line = -1.7, cex = 0.75, adj = 0.98
    )
    graphics::legend(
      "bottomright",
      legend = c("PDP", "CSPD"),
      col = c("grey20", "#0072B2"),
      lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8
    )
  }
  path
}

plot_s1_support <- function(simulation, path = file.path("results", "s1-support-diagnostics.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  rho_values <- sort(unique(simulation$rho))
  grDevices::png(path, width = 1800, height = 1200, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  for (rho in rho_values) {
    subset <- simulation[simulation$rho == rho, ]
    graphics::plot(
      subset$z, subset$isc,
      type = "l", ylim = c(0, 1), lwd = 2, col = "#009E73",
      xlab = expression(z), ylab = "Coverage"
    )
    graphics::lines(subset$z, subset$per, lwd = 2, col = "#D55E00", lty = 2)
    graphics::abline(h = unique(subset$gamma), col = "#0072B2", lwd = 2, lty = 3)
    graphics::legend(
      "bottom",
      legend = c("ISC", "PER", expression(gamma(G))),
      col = c("#009E73", "#D55E00", "#0072B2"),
      lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.75
    )
  }
  path
}
