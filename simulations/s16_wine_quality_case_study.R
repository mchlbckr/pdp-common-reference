# Case study S16: held-out red-wine quality prediction (UCI Wine Quality data).

load_s16_wine <- function(path = file.path("data", "winequality-red.csv")) {
  candidates <- c(path, file.path("..", "..", path))
  available <- candidates[file.exists(candidates)]
  if (length(available) == 0L) stop("wine-quality data not found at ", path, call. = FALSE)
  utils::read.csv(available[[1]], sep = ";", check.names = TRUE)
}

s16_fit_models <- function(training) {
  additive <- mgcv::gam(
    quality ~ s(alcohol, k = 8) + s(volatile.acidity, k = 8) +
      s(sulphates, k = 8) + s(density, k = 8) +
      s(total.sulfur.dioxide, k = 8) + s(chlorides, k = 8),
    data = training, method = "REML"
  )
  interaction <- mgcv::gam(
    quality ~ s(alcohol, k = 8) + s(volatile.acidity, k = 8) +
      s(sulphates, k = 8) + s(density, k = 8) +
      s(total.sulfur.dioxide, k = 8) + s(chlorides, k = 8) +
      ti(alcohol, density, k = c(8, 8)),
    data = training, method = "REML"
  )
  boosted <- gbm::gbm(
    quality ~ alcohol + volatile.acidity + sulphates + density +
      total.sulfur.dioxide + chlorides,
    data = training, distribution = "gaussian", n.trees = 300L,
    interaction.depth = 3L, shrinkage = 0.03, n.minobsinnode = 15L,
    verbose = FALSE
  )
  list(additive_gam = additive, interaction_gam = interaction, boosted_tree = boosted)
}

s16_predictor <- function(model) {
  if (inherits(model, "gbm")) {
    return(function(newdata) as.numeric(stats::predict(
      model, newdata = newdata, n.trees = model$n.trees
    )))
  }
  function(newdata) as.numeric(stats::predict(model, newdata = newdata))
}

run_s16 <- function(alpha = 0.10, relaxed_epsilon = 0.10, grid_points = 41L, seed = 20260922L) {
  set.seed(seed)
  data <- load_s16_wine()
  index <- sample.int(nrow(data))
  n_training <- floor(0.60 * nrow(data))
  n_calibration <- floor(0.20 * nrow(data))
  training <- data[index[seq_len(n_training)], ]
  calibration <- data[index[n_training + seq_len(n_calibration)], ]
  evaluation <- data[index[(n_training + n_calibration + 1L):nrow(data)], ]
  feature_names <- setdiff(names(data), "quality")
  models <- s16_fit_models(training)
  validity_fit <- fit_locally_scaled_validity(
    training[, feature_names], calibration[, feature_names], "alcohol", alpha = alpha
  )
  is_valid <- function(newdata) predict_locally_scaled_validity(validity_fit, newdata)
  # The central 30th--70th percentile grid is declared before evaluation.  A
  # wider grid has zero common mass for this validity criterion,
  # which is itself reported in the manuscript as a grid-sensitivity finding.
  grid <- seq(stats::quantile(training$alcohol, 0.30),
              stats::quantile(training$alcohol, 0.70), length.out = grid_points)
  evaluation_features <- evaluation[, feature_names]
  isc <- estimate_isc(evaluation_features, "alcohol", grid, is_valid)
  evaluate_model <- function(model_name) {
    predict_function <- s16_predictor(models[[model_name]])
    pdp <- estimate_pdp(evaluation_features, "alcohol", grid, predict_function)
    cspd <- estimate_cspd(evaluation_features, "alcohol", grid, is_valid, predict_function)
    relaxed <- estimate_relaxed_cspd(
      evaluation_features, "alcohol", grid, is_valid, predict_function,
      epsilon = relaxed_epsilon, tolerance = 2e-5
    )
    data.frame(
      z = grid, pdp = pdp$pdp, cspd = cspd$cspd, isc = isc$isc,
      relaxed_cspd = relaxed$relaxed_cspd, relaxed_isc = relaxed$weighted_coverage,
      gamma = cspd$gamma, n_common = cspd$n_common,
      relaxed_epsilon = relaxed_epsilon, relaxed_min_coverage = relaxed$min_coverage,
      relaxed_mean_coverage = relaxed$mean_coverage,
      relaxed_simultaneous_valid_mass = relaxed$simultaneous_valid_mass,
      relaxed_kl = relaxed$kl, relaxed_ess = relaxed$effective_sample_size,
      relaxed_max_weight = relaxed$max_weight,
      active_constraints = relaxed$active_constraints,
      active_affine_dimension = relaxed$active_affine_dimension,
      active_licq = relaxed$active_licq,
      binding_constraints = relaxed$binding_constraints,
      strict_complementarity = relaxed$strict_complementarity,
      observed_patterns = relaxed$observed_patterns,
      median_pattern_size = relaxed$median_pattern_size,
      min_pattern_size = relaxed$min_pattern_size,
      singleton_pattern_fraction = relaxed$singleton_pattern_fraction,
      epsilon_min = relaxed$epsilon_min,
      weight_prediction_correlation = relaxed$weight_prediction_correlation,
      max_abs_weight_prediction_correlation = max(abs(relaxed$weight_prediction_correlation)),
      reference_covariance = relaxed$reference_covariance,
      covariance_identity_error = relaxed$covariance_identity_error,
      max_covariance_identity_error = max(relaxed$covariance_identity_error),
      covariance_path_range = diff(range(relaxed$reference_covariance)),
      centered_hard_pdp_distance = centered_curve_distance(cspd$cspd, pdp$pdp),
      centered_relaxed_pdp_distance = centered_curve_distance(relaxed$relaxed_cspd, pdp$pdp),
      centered_relaxed_hard_distance = centered_curve_distance(relaxed$relaxed_cspd, cspd$cspd),
      heldout_mae = mean(abs(evaluation$quality - predict_function(evaluation_features))),
      n_training = nrow(training), n_calibration = nrow(calibration),
      n_evaluation = nrow(evaluation), model_name = model_name,
      stringsAsFactors = FALSE
    )
  }
  evaluated <- lapply(names(models), evaluate_model)
  names(evaluated) <- names(models)
  candidate_summary <- do.call(rbind, lapply(evaluated, function(result) {
    result[1, c("model_name", "heldout_mae", "centered_hard_pdp_distance",
                "centered_relaxed_pdp_distance", "centered_relaxed_hard_distance")]
  }))
  # The flexible comparator is declared by model class, not selected for a
  # large displayed RCPD--PD difference on the evaluation set.
  selected_name <- "boosted_tree"
  keep <- c("additive_gam", selected_name)
  output <- do.call(rbind, evaluated[keep])
  output$model_case <- ifelse(
    output$model_name == "additive_gam", "Additive GAM",
    ifelse(output$model_name == "interaction_gam", "GAM with alcohol-density interaction",
           "Gradient boosting")
  )
  widths <- seq(0.2, 0.9, by = 0.1)
  profile_alphas <- c(0.01, 0.05, 0.10, 0.20)
  attr(output, "support_profile") <- do.call(rbind, lapply(profile_alphas, function(profile_alpha) {
    profile_fit <- fit_locally_scaled_validity(
      training[, feature_names], calibration[, feature_names], "alcohol", alpha = profile_alpha
    )
    profile_valid <- function(newdata) predict_locally_scaled_validity(profile_fit, newdata)
    do.call(rbind, lapply(widths, function(width) {
      profile_grid <- seq(stats::quantile(training$alcohol, (1 - width) / 2),
                          stats::quantile(training$alcohol, 1 - (1 - width) / 2),
                          length.out = grid_points)
      data.frame(alpha = profile_alpha, central_fraction = width,
                 gamma = mean(common_support_indicator(
                   evaluation_features, "alcohol", profile_grid, profile_valid
                 )))
    }))
  }))
  attr(output, "sample_focal") <- evaluation_features$alcohol
  attr(output, "evaluation_contexts") <- evaluation_features
  attr(output, "candidate_summary") <- candidate_summary
  attr(output, "selected_nonadditive_model") <- selected_name
  attr(output, "evaluation_validity") <- validity_matrix(
    evaluation_features, "alcohol", grid, is_valid
  )
  all_prediction_matrices <- lapply(names(models), function(model_name) {
    prediction_matrix(
      evaluation_features, "alcohol", grid,
      s16_predictor(models[[model_name]])
    )
  })
  names(all_prediction_matrices) <- names(models)
  attr(output, "all_prediction_matrices") <- all_prediction_matrices
  attr(output, "prediction_matrices") <-
    all_prediction_matrices[c("additive_gam", selected_name)]
  output
}

plot_s16 <- function(case_study, path = file.path("results", "s16-wine-quality-case-study.png")) {
  curve_colors <- c("PDP" = pdp_palette[["ink"]], "Hard CSPD" = pdp_palette[["orange"]],
                    "RCPD (epsilon = 0.10)" = pdp_palette[["blue"]])
  support_colors <- c("Original reference P" = pdp_palette[["green"]],
                      "Relaxed reference Q* (epsilon = 0.10)" = pdp_palette[["blue"]])
  profile <- attr(case_study, "support_profile")
  rug <- data.frame(z = attr(case_study, "sample_focal"))
  plots <- list()
  for (model_label in unique(case_study$model_case)) {
    case <- case_study[case_study$model_case == model_label, , drop = FALSE]
    curves <- rbind(
      data.frame(z = case$z, value = case$pdp - mean(case$pdp), curve = "PDP"),
      data.frame(z = case$z, value = case$cspd - mean(case$cspd), curve = "Hard CSPD"),
      data.frame(z = case$z, value = case$relaxed_cspd - mean(case$relaxed_cspd),
                 curve = "RCPD (epsilon = 0.10)")
    )
    curves$curve <- factor(curves$curve,
                           levels = c("PDP", "Hard CSPD", "RCPD (epsilon = 0.10)"))
    visible_rug <- rug[rug$z >= min(curves$z) & rug$z <= max(curves$z), , drop = FALSE]
    p_curves <- ggplot2::ggplot(
      curves, ggplot2::aes(z, value, color = curve, linetype = curve)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_rug(data = visible_rug, ggplot2::aes(x = z), inherit.aes = FALSE,
                        sides = "b", alpha = 0.18, color = pdp_palette[["mid_grey"]]) +
      ggplot2::annotate(
        "text", x = Inf, y = Inf,
        label = sprintf("MAE %.3f; RCPD-PD shape %.3f",
                        unique(case$heldout_mae),
                        unique(case$centered_relaxed_pdp_distance)),
        hjust = 1.05, vjust = 1.3, size = 2.6
      ) +
      ggplot2::scale_color_manual(values = curve_colors) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
      ggplot2::labs(title = model_label, x = "Alcohol (% by volume)",
                    y = "Centered predicted quality", color = NULL, linetype = NULL) +
      theme_pdp(base_size = 9)
    support <- rbind(
      data.frame(z = case$z, value = case$isc, diagnostic = "Original reference P"),
      data.frame(z = case$z, value = case$relaxed_isc,
                 diagnostic = "Relaxed reference Q* (epsilon = 0.10)")
    )
    support$diagnostic <- factor(
      support$diagnostic,
      levels = c("Original reference P", "Relaxed reference Q* (epsilon = 0.10)")
    )
    p_support <- ggplot2::ggplot(
      support, ggplot2::aes(z, value, color = diagnostic, linetype = diagnostic)
    ) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::scale_color_manual(values = support_colors) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = "Alcohol (% by volume)", y = "Pointwise support coverage",
                    color = NULL, linetype = NULL) +
      theme_pdp(base_size = 9)
    p_profile <- ggplot2::ggplot(
      profile, ggplot2::aes(central_fraction, gamma, color = factor(alpha),
                            linetype = factor(alpha), group = alpha)
    ) +
      ggplot2::geom_hline(yintercept = 0.1, linetype = "dotted",
                          color = pdp_palette[["mid_grey"]]) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.7) +
      ggplot2::scale_color_manual(
        values = c("0.01" = pdp_palette[["green"]], "0.05" = pdp_palette[["orange"]],
                   "0.1" = pdp_palette[["blue"]], "0.2" = pdp_palette[["purple"]])
      ) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = "Central grid fraction", y = expression(hat(gamma)(G)),
                    color = expression(alpha), linetype = expression(alpha)) +
      theme_pdp(base_size = 9)
    plots <- c(plots, list(p_curves, p_support, p_profile))
  }
  save_pdp_panel_set(plots, path, widths = 4.2, heights = 3.7)
}
