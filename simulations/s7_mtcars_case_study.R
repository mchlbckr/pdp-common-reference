# Case study S7: a reproducible observational-data illustration using mtcars.

run_s7 <- function(k = 4L, alpha = 0.10, grid_points = 41L) {
  raw <- datasets::mtcars
  model_data <- raw[, c("mpg", "wt", "hp", "qsec")]
  reference_data <- model_data[, c("wt", "hp", "qsec")]
  fitted_model <- stats::lm(mpg ~ wt * hp + qsec, data = model_data)
  predict_function <- function(newdata) stats::predict(fitted_model, newdata = newdata)
  validity_fit <- fit_knn_validity(reference_data, k = k, alpha = alpha)
  is_valid <- function(newdata) predict_knn_validity(validity_fit, newdata)

  grid <- seq(
    stats::quantile(reference_data$wt, 0.05),
    stats::quantile(reference_data$wt, 0.95),
    length.out = grid_points
  )
  pdp <- estimate_pdp(reference_data, "wt", grid, predict_function)
  cspd <- estimate_cspd(reference_data, "wt", grid, is_valid, predict_function)
  isc <- estimate_isc(reference_data, "wt", grid, is_valid)

  half_width <- seq(0.05, 0.45, by = 0.05)
  grid_sensitivity <- do.call(rbind, lapply(half_width, function(width) {
    candidate_grid <- seq(
      stats::quantile(reference_data$wt, width),
      stats::quantile(reference_data$wt, 1 - width),
      length.out = grid_points
    )
    candidate_cspd <- estimate_cspd(reference_data, "wt", candidate_grid, is_valid, predict_function)
    data.frame(
      lower_quantile = width,
      upper_quantile = 1 - width,
      grid_width = 1 - 2 * width,
      gamma = unique(candidate_cspd$gamma),
      n_common = unique(candidate_cspd$n_common)
    )
  }))

  list(
    curves = cbind(pdp, cspd = cspd$cspd, isc = isc$isc, per = isc$per),
    grid_sensitivity = grid_sensitivity,
    summary = data.frame(
      n = nrow(reference_data), k = k, alpha = alpha,
      gamma = unique(cspd$gamma), n_common = unique(cspd$n_common),
      minimum_isc = min(isc$isc), maximum_per = max(isc$per)
    )
  )
}

plot_s7 <- function(case_study, path = file.path("results", "s7-mtcars-case-study.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  curves <- case_study$curves
  sensitivity <- case_study$grid_sensitivity
  sensitivity <- sensitivity[order(sensitivity$grid_width), , drop = FALSE]
  summary <- case_study$summary
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

  graphics::plot(curves$z, curves$pdp, type = "l", lwd = 2, col = "grey20",
                 xlab = "Vehicle weight (1000 lb)", ylab = "Predicted mpg")
  graphics::lines(curves$z, curves$cspd, lwd = 2, lty = 2, col = "#0072B2")
  graphics::legend("topright", legend = c("PDP", "CSPD"), col = c("grey20", "#0072B2"),
                   lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)

  graphics::plot(curves$z, curves$isc, type = "l", ylim = c(0, 1), lwd = 2,
                 col = "#009E73", xlab = "Vehicle weight (1000 lb)",
                 ylab = "Support coverage")
  graphics::lines(curves$z, curves$per, lwd = 2, lty = 2, col = "#D55E00")
  graphics::abline(h = summary$gamma, lwd = 2, lty = 3, col = "#0072B2")
  graphics::legend("bottom", legend = c("ISC", "PER", expression(gamma(G))),
                   col = c("#009E73", "#D55E00", "#0072B2"), lty = c(1, 2, 3),
                   lwd = 2, bty = "n", cex = 0.8)

  graphics::plot(sensitivity$grid_width, sensitivity$gamma, type = "b", pch = 19,
                 lwd = 2, ylim = c(0, 1), col = "#0072B2", xlab = "Central quantile-grid width",
                 ylab = expression(gamma(G)))
  graphics::text(sensitivity$grid_width, sensitivity$gamma, labels = sensitivity$n_common,
                 pos = 3, cex = 0.7, col = "grey25")
  graphics::legend("bottomleft", legend = "Point labels: common contexts", bty = "n", cex = 0.7)
  path
}
