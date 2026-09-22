# Simulation S12: support-learner robustness on a nonlinear manifold.

run_s12 <- function(
  n_training = 150L,
  n_calibration = 150L,
  n_evaluation = 300L,
  repetitions = 40L,
  sigma = 0.12,
  alpha = 0.05,
  k = 15L,
  grid = seq(-0.4, 0.4, length.out = 15L),
  seed = 20260920L
) {
  grid <- assert_grid(grid)
  oracle_validity <- function(newdata) {
    abs(newdata$x2 - newdata$x1^2) <= stats::qnorm(1 - alpha / 2) * sigma
  }
  set.seed(seed)
  raw <- do.call(rbind, lapply(seq_len(repetitions), function(replication) {
    training <- simulate_s2_data(n_training, sigma)
    calibration <- simulate_s2_data(n_calibration, sigma)
    evaluation <- simulate_s2_data(n_evaluation, sigma)
    linear_fit <- fit_split_conformal_validity(training, calibration, "x1", alpha)
    knn_fit <- fit_knn_validity(training, k = k, alpha = alpha,
                                calibration_data = calibration)
    learners <- list(
      "linear split-conformal" = function(newdata) predict_validity(linear_fit, newdata),
      "local kNN" = function(newdata) predict_knn_validity(knn_fit, newdata)
    )
    oracle_isc <- estimate_isc(evaluation, "x1", grid, oracle_validity)
    oracle_cspd <- estimate_cspd(evaluation, "x1", grid, oracle_validity, s2_prediction)
    do.call(rbind, lapply(names(learners), function(learner) {
      validity <- learners[[learner]]
      isc <- estimate_isc(evaluation, "x1", grid, validity)
      cspd <- estimate_cspd(evaluation, "x1", grid, validity, s2_prediction)
      data.frame(
        learner = learner,
        isc_mae = mean(abs(isc$isc - oracle_isc$isc)),
        cspd_mae = mean(abs(cspd$cspd - oracle_cspd$cspd), na.rm = TRUE),
        gamma_error = abs(unique(cspd$gamma) - unique(oracle_cspd$gamma)),
        gamma = unique(cspd$gamma),
        oracle_gamma = unique(oracle_cspd$gamma)
      )
    }))
  }))
  split_raw <- split(raw, raw$learner)
  do.call(rbind, lapply(split_raw, function(value) {
    data.frame(
      learner = unique(value$learner),
      isc_mae = mean(value$isc_mae),
      cspd_mae = mean(value$cspd_mae),
      gamma_error = mean(value$gamma_error),
      gamma = mean(value$gamma),
      oracle_gamma = mean(value$oracle_gamma),
      stringsAsFactors = FALSE
    )
  }))
}

plot_s12 <- function(simulation, path = file.path("results", "s12-validity-learner-robustness.png")) {
  display <- c("linear split-conformal" = "Linear split-\nconformal",
               "local kNN" = "Local kNN")
  colors <- c("Linear split-\nconformal" = pdp_palette[["orange"]],
              "Local kNN" = pdp_palette[["blue"]])
  make_panel <- function(measure, y_label) {
    data <- data.frame(
      learner = factor(unname(display[simulation$learner]), levels = unname(display)),
      value = simulation[[measure]]
    )
    ggplot2::ggplot(data, ggplot2::aes(learner, value, fill = learner)) +
      ggplot2::geom_col(width = 0.64, color = pdp_palette[["ink"]], linewidth = 0.35) +
      ggplot2::scale_fill_manual(values = colors, guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
      ggplot2::labs(x = NULL, y = y_label) +
      theme_pdp() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                     legend.position = "none")
  }
  panels <- list(
    make_panel("isc_mae", "Mean absolute ISC error"),
    make_panel("cspd_mae", "Mean absolute CSPD error"),
    make_panel("gamma_error", expression("Mean absolute error in " * gamma(G)))
  )
  save_pdp_panel_set(panels, path, widths = 4.2, heights = 4.0)
}
