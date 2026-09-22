# Case study S15: Boston housing benchmark (MASS), excluding the race-derived field.

load_s15_boston <- function() {
  data_environment <- new.env(parent = emptyenv())
  utils::data("Boston", package = "MASS", envir = data_environment)
  data_environment$Boston[, c("medv", "lstat", "rm", "ptratio", "nox", "dis")]
}

run_s15 <- function(k = 25L, alpha = 0.10, grid_points = 41L) {
  model_data <- load_s15_boston()
  reference_data <- model_data[, setdiff(names(model_data), "medv")]
  fitted_model <- rpart::rpart(medv ~ ., data = model_data,
                               control = rpart::rpart.control(cp = 0.01, minsplit = 20L))
  predict_function <- function(newdata) as.numeric(stats::predict(fitted_model, newdata = newdata))
  validity_fit <- fit_knn_validity(reference_data, k = k, alpha = alpha)
  is_valid <- function(newdata) predict_knn_validity(validity_fit, newdata)
  grid <- seq(stats::quantile(reference_data$lstat, 0.05),
              stats::quantile(reference_data$lstat, 0.95), length.out = grid_points)
  pdp <- estimate_pdp(reference_data, "lstat", grid, predict_function)
  cspd <- estimate_cspd(reference_data, "lstat", grid, is_valid, predict_function)
  isc <- estimate_isc(reference_data, "lstat", grid, is_valid)
  data.frame(z = grid, pdp = pdp$pdp, cspd = cspd$cspd, isc = isc$isc, per = isc$per,
             gamma = cspd$gamma, n_common = cspd$n_common)
}

plot_s15 <- function(case_study, path = file.path("results", "s15-boston-case-study.png")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 1800, height = 900, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old_par); grDevices::dev.off() }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0))
  y_limits <- range(c(case_study$pdp, case_study$cspd), na.rm = TRUE)
  y_padding <- 0.04 * diff(y_limits)
  graphics::plot(case_study$z, case_study$pdp, type = "l", lwd = 2, col = "grey20",
                 ylim = y_limits + c(-y_padding, y_padding),
                 xlab = "Lower-status percentage", ylab = "Predicted median home value")
  graphics::lines(case_study$z, case_study$cspd, lwd = 2, lty = 2, col = "#0072B2")
  graphics::legend("topright", legend = c("PDP", "CSPD"), col = c("grey20", "#0072B2"),
                   lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  graphics::plot(case_study$z, case_study$isc, type = "l", ylim = c(0, 1), lwd = 2,
                 col = "#009E73", xlab = "Lower-status percentage", ylab = "Support coverage")
  graphics::lines(case_study$z, case_study$per, lwd = 2, lty = 2, col = "#D55E00")
  graphics::abline(h = unique(case_study$gamma), lwd = 2, lty = 3, col = "#0072B2")
  graphics::legend("bottom", legend = c("ISC", "PER", expression(gamma(G))),
                   col = c("#009E73", "#D55E00", "#0072B2"), lty = c(1, 2, 3),
                   lwd = 2, bty = "n", cex = 0.8)
  path
}
