# Case study S13: airquality with a nonlinear fitted ozone model.

run_s13 <- function(k = 10L, alpha = 0.10, grid_points = 41L) {
  model_data <- stats::na.omit(datasets::airquality[, c("Ozone", "Temp", "Wind", "Solar.R")])
  reference_data <- model_data[, c("Temp", "Wind", "Solar.R")]
  fitted_model <- stats::lm(Ozone ~ stats::poly(Temp, 2) + Wind * Temp + Solar.R, data = model_data)
  predict_function <- function(newdata) stats::predict(fitted_model, newdata = newdata)
  validity_fit <- fit_knn_validity(reference_data, k = k, alpha = alpha)
  is_valid <- function(newdata) predict_knn_validity(validity_fit, newdata)
  grid <- seq(stats::quantile(reference_data$Temp, 0.05),
              stats::quantile(reference_data$Temp, 0.95), length.out = grid_points)
  pdp <- estimate_pdp(reference_data, "Temp", grid, predict_function)
  cspd <- estimate_cspd(reference_data, "Temp", grid, is_valid, predict_function)
  isc <- estimate_isc(reference_data, "Temp", grid, is_valid)
  data.frame(
    z = grid, pdp = pdp$pdp, cspd = cspd$cspd, isc = isc$isc, per = isc$per,
    gamma = cspd$gamma, n_common = cspd$n_common
  )
}

plot_s13 <- function(case_study, path = file.path("results", "s13-airquality-case-study.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0))
  graphics::plot(case_study$z, case_study$pdp, type = "l", lwd = 2, col = "grey20",
                 xlab = "Temperature (F)", ylab = "Predicted ozone")
  graphics::lines(case_study$z, case_study$cspd, lwd = 2, lty = 2, col = "#0072B2")
  graphics::legend("topleft", legend = c("PDP", "CSPD"), col = c("grey20", "#0072B2"),
                   lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  graphics::plot(case_study$z, case_study$isc, type = "l", ylim = c(0, 1), lwd = 2,
                 col = "#009E73", xlab = "Temperature (F)", ylab = "Support coverage")
  graphics::lines(case_study$z, case_study$per, lwd = 2, lty = 2, col = "#D55E00")
  graphics::abline(h = unique(case_study$gamma), lwd = 2, lty = 3, col = "#0072B2")
  graphics::legend("bottomleft", inset = c(0.02, 0.02), legend = c("ISC", "PER", expression(gamma(G))),
                   col = c("#009E73", "#D55E00", "#0072B2"), lty = c(1, 2, 3),
                   lwd = 2, bty = "n", cex = 0.8)
  path
}
