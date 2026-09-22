# Held-out, mixed-type tabular applications using fixed one-hot feature coding.

s17_encode_predictors <- function(data, outcome, drop = character()) {
  predictors <- data[, setdiff(names(data), c(outcome, drop)), drop = FALSE]
  character_columns <- vapply(predictors, is.character, logical(1))
  predictors[character_columns] <- lapply(predictors[character_columns], factor)
  logical_columns <- vapply(predictors, is.logical, logical(1))
  predictors[logical_columns] <- lapply(predictors[logical_columns], as.integer)
  design <- stats::model.matrix(~ . - 1, data = predictors, na.action = stats::na.pass)
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  as.data.frame(design, check.names = TRUE)
}

s17_data_path <- function(...) {
  relative <- file.path(...)
  candidates <- c(relative, file.path("..", "..", relative))
  available <- candidates[file.exists(candidates)]
  if (length(available) == 0L) stop("data file not found: ", relative, call. = FALSE)
  available[[1]]
}

s17_split <- function(data, seed) {
  set.seed(seed)
  index <- sample.int(nrow(data))
  n_training <- floor(0.60 * nrow(data))
  n_calibration <- floor(0.20 * nrow(data))
  list(training = index[seq_len(n_training)],
       calibration = index[n_training + seq_len(n_calibration)],
       evaluation = index[(n_training + n_calibration + 1L):nrow(data)])
}

s17_summarize_curves <- function(training, calibration, evaluation, focal, prediction,
                                 alpha = 0.10, relaxed_epsilon = 0.10,
                                 grid_quantiles = c(0.30, 0.70)) {
  validity_fit <- suppressWarnings(fit_locally_scaled_validity(training, calibration, focal, alpha = alpha))
  is_valid <- function(newdata) suppressWarnings(predict_locally_scaled_validity(validity_fit, newdata))
  grid <- seq(stats::quantile(training[[focal]], grid_quantiles[[1]]),
              stats::quantile(training[[focal]], grid_quantiles[[2]]), length.out = 41L)
  pdp <- estimate_pdp(evaluation, focal, grid, prediction)
  cspd <- estimate_cspd(evaluation, focal, grid, is_valid, prediction)
  relaxed <- estimate_relaxed_cspd(evaluation, focal, grid, is_valid, prediction,
                                   epsilon = relaxed_epsilon, tolerance = 2e-5)
  isc <- estimate_isc(evaluation, focal, grid, is_valid)
  output <- data.frame(z = grid, pdp = pdp$pdp, cspd = cspd$cspd,
             relaxed_cspd = relaxed$relaxed_cspd, relaxed_isc = relaxed$weighted_coverage,
             isc = isc$isc, per = isc$per,
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
             centered_relaxed_pdp_distance = centered_curve_distance(relaxed$relaxed_cspd, pdp$pdp))
  widths <- seq(0.2, 0.9, by = 0.1)
  profile_alphas <- c(0.01, 0.05, 0.10, 0.20)
  attr(output, "support_profile") <- do.call(rbind, lapply(profile_alphas, function(profile_alpha) {
    profile_fit <- suppressWarnings(fit_locally_scaled_validity(
      training, calibration, focal, alpha = profile_alpha
    ))
    profile_valid <- function(newdata) {
      suppressWarnings(predict_locally_scaled_validity(profile_fit, newdata))
    }
    do.call(rbind, lapply(widths, function(width) {
      profile_grid <- seq(stats::quantile(training[[focal]], (1 - width) / 2),
                          stats::quantile(training[[focal]], 1 - (1 - width) / 2),
                          length.out = 41L)
      data.frame(alpha = profile_alpha, central_fraction = width,
                 gamma = mean(common_support_indicator(
                   evaluation, focal, profile_grid, profile_valid
                 )))
    }))
  }))
  attr(output, "sample_focal") <- evaluation[[focal]]
  attr(output, "evaluation_validity") <- validity_matrix(
    evaluation, focal, grid, is_valid
  )
  attr(output, "prediction_matrix") <- prediction_matrix(
    evaluation, focal, grid, prediction
  )
  attr(output, "observed_validity") <- is_valid(evaluation)
  attr(output, "observed_predictions") <- prediction(evaluation)
  output
}

run_s17_bank <- function(seed = 20261120L) {
  raw <- utils::read.csv(s17_data_path("data", "uci", "bank-extract", "bank-full.csv"), sep = ";")
  # Call duration is unavailable before contact and is excluded to avoid target leakage.
  encoded <- s17_encode_predictors(raw, outcome = "y", drop = "duration")
  split <- s17_split(encoded, seed)
  training <- encoded[split$training, , drop = FALSE]
  calibration <- encoded[split$calibration, , drop = FALSE]
  evaluation <- encoded[split$evaluation, , drop = FALSE]
  response <- raw$y == "yes"
  model <- randomForest::randomForest(x = training, y = factor(response[split$training]),
                                      ntree = 100L, mtry = max(1L, floor(sqrt(ncol(training)))))
  prediction <- function(newdata) as.numeric(stats::predict(model, newdata = newdata, type = "prob")[, "TRUE"])
  output <- s17_summarize_curves(training, calibration, evaluation, "age", prediction)
  support_profile <- attr(output, "support_profile")
  sample_focal <- attr(output, "sample_focal")
  evaluation_validity <- attr(output, "evaluation_validity")
  evaluation_predictions <- attr(output, "prediction_matrix")
  observed_validity <- attr(output, "observed_validity")
  observed_predictions <- attr(output, "observed_predictions")
  joint_fit <- fit_joint_gaussian_validity(training, calibration, alpha = 0.10)
  joint_valid <- function(newdata) predict_joint_gaussian_validity(joint_fit, newdata)
  joint_gamma <- mean(common_support_indicator(evaluation, "age", output$z, joint_valid))
  joint_mean_isc <- mean(validity_matrix(evaluation, "age", output$z, joint_valid))
  output <- transform(output, dataset = "Bank marketing", model = "Random forest", focal = "Age",
            conditional_gamma = gamma, joint_gamma = joint_gamma,
            joint_mean_isc = joint_mean_isc,
            heldout_brier = mean((response[split$evaluation] - prediction(evaluation))^2),
            n_training = nrow(training), n_calibration = nrow(calibration), n_evaluation = nrow(evaluation))
  attr(output, "support_profile") <- support_profile
  attr(output, "sample_focal") <- sample_focal
  attr(output, "evaluation_validity") <- evaluation_validity
  attr(output, "prediction_matrix") <- evaluation_predictions
  attr(output, "observed_validity") <- observed_validity
  attr(output, "observed_predictions") <- observed_predictions
  attr(output, "evaluation_response") <- response[split$evaluation]
  attr(output, "score_comparison") <- data.frame(
    score = c("Conditional locally scaled", "Joint Gaussian"),
    mean_isc = c(mean(output$isc), joint_mean_isc),
    gamma = c(unique(output$gamma), joint_gamma), stringsAsFactors = FALSE
  )
  output
}

run_s17_shoppers <- function(seed = 20261121L) {
  raw <- utils::read.csv(s17_data_path("data", "uci", "online_shoppers_intention.csv"))
  encoded <- s17_encode_predictors(raw, outcome = "Revenue")
  split <- s17_split(encoded, seed)
  training <- encoded[split$training, , drop = FALSE]
  calibration <- encoded[split$calibration, , drop = FALSE]
  evaluation <- encoded[split$evaluation, , drop = FALSE]
  response <- as.integer(raw$Revenue)
  training_model <- cbind(training, outcome = response[split$training])
  model <- gbm::gbm(outcome ~ ., data = training_model, distribution = "bernoulli", n.trees = 150L,
                    interaction.depth = 2L, shrinkage = 0.03, n.minobsinnode = 20L, verbose = FALSE)
  prediction <- function(newdata) stats::plogis(as.numeric(stats::predict(model, newdata = newdata, n.trees = model$n.trees)))
  output <- s17_summarize_curves(training, calibration, evaluation, "ProductRelated_Duration", prediction)
  support_profile <- attr(output, "support_profile")
  sample_focal <- attr(output, "sample_focal")
  evaluation_validity <- attr(output, "evaluation_validity")
  evaluation_predictions <- attr(output, "prediction_matrix")
  observed_validity <- attr(output, "observed_validity")
  observed_predictions <- attr(output, "observed_predictions")
  output <- transform(output, dataset = "Online shoppers", model = "Gradient boosting", focal = "Product-related duration",
            heldout_brier = mean((response[split$evaluation] - prediction(evaluation))^2),
            n_training = nrow(training), n_calibration = nrow(calibration), n_evaluation = nrow(evaluation))
  attr(output, "support_profile") <- support_profile
  attr(output, "sample_focal") <- sample_focal
  attr(output, "evaluation_validity") <- evaluation_validity
  attr(output, "prediction_matrix") <- evaluation_predictions
  attr(output, "observed_validity") <- observed_validity
  attr(output, "observed_predictions") <- observed_predictions
  attr(output, "evaluation_response") <- response[split$evaluation]
  output
}

plot_s17 <- function(case_studies, path = file.path("results", "s17-large-mixed-case-studies.png")) {
  plots <- list()
  for (case in case_studies) {
    curves <- rbind(
      data.frame(z = case$z, value = case$pdp - mean(case$pdp), curve = "PDP"),
      data.frame(z = case$z, value = case$cspd - mean(case$cspd), curve = "Hard CSPD"),
      data.frame(z = case$z, value = case$relaxed_cspd - mean(case$relaxed_cspd), curve = "RCPD (epsilon = 0.10)")
    )
    curves$curve <- factor(curves$curve,
                           levels = c("PDP", "Hard CSPD", "RCPD (epsilon = 0.10)"))
    rug <- data.frame(z = attr(case, "sample_focal"))
    rug <- rug[rug$z >= min(curves$z) & rug$z <= max(curves$z), , drop = FALSE]
    curve_colors <- c("PDP" = pdp_palette[["ink"]],
                      "Hard CSPD" = pdp_palette[["orange"]],
                      "RCPD (epsilon = 0.10)" = pdp_palette[["blue"]])
    p_curves <- ggplot2::ggplot(curves, ggplot2::aes(z, value, color = curve,
                                                     linetype = curve)) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::geom_rug(data = rug,
                        ggplot2::aes(x = z), inherit.aes = FALSE, sides = "b",
                        alpha = 0.14, color = pdp_palette[["mid_grey"]]) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = sprintf("Held-out Brier: %.3f", unique(case$heldout_brier)),
                        hjust = 1.05, vjust = 1.3, size = 2.6) +
      ggplot2::scale_color_manual(values = curve_colors) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
      ggplot2::labs(title = unique(case$dataset), x = unique(case$focal),
                    y = "Centered predicted probability", color = NULL, linetype = NULL) +
      theme_pdp(base_size = 9)
    support <- rbind(
      data.frame(z = case$z, value = case$isc,
                 diagnostic = "Original reference P"),
      data.frame(z = case$z, value = case$relaxed_isc,
                 diagnostic = "Relaxed reference Q* (epsilon = 0.10)")
    )
    support$diagnostic <- factor(support$diagnostic,
                                 levels = c("Original reference P", "Relaxed reference Q* (epsilon = 0.10)"))
    support_colors <- c("Original reference P" = pdp_palette[["green"]],
                        "Relaxed reference Q* (epsilon = 0.10)" = pdp_palette[["blue"]])
    p_support <- ggplot2::ggplot(support, ggplot2::aes(z, value, color = diagnostic,
                                                       linetype = diagnostic)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::scale_color_manual(values = support_colors) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = unique(case$focal), y = "Pointwise support coverage",
                    color = NULL, linetype = NULL) +
      theme_pdp(base_size = 9)
    profile <- attr(case, "support_profile")
    p_profile <- ggplot2::ggplot(
      profile, ggplot2::aes(central_fraction, gamma, color = factor(alpha),
                            linetype = factor(alpha), group = alpha)
    ) +
      ggplot2::geom_hline(yintercept = 0.1, linetype = "dotted",
                          color = pdp_palette[["mid_grey"]]) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::scale_color_manual(
        values = c("0.01" = pdp_palette[["green"]], "0.05" = pdp_palette[["orange"]],
                   "0.1" = pdp_palette[["blue"]], "0.2" = pdp_palette[["purple"]])
      ) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = "Central grid fraction", y = expression(hat(gamma)(G)),
                    color = expression(alpha), linetype = expression(alpha)) +
      theme_pdp(base_size = 9)
    plots <- c(plots, list(p_curves, p_support, p_profile))
  }
  save_pdp_panel_set(plots, path, widths = 4.2, heights = 3.7)
}
