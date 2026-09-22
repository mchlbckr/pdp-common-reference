# Confirmatory model-class check for the California-housing Longitude effect.
#
# This script is deliberately separate from the paper pipeline. Its design is
# fixed before execution: 30 fresh splits, random forest and gradient boosting,
# two validity learners, epsilon in {0.10, 0.15, 0.20}, alpha = 0.10, the
# central 20% of training-set Longitude, and K = 21. The primary gate requires
# at least one learner/epsilon setting to pass all criteria for both model
# classes; no manuscript result is changed automatically. Boosting capacity was
# fixed on the excluded smoke-test seed 20269991: 150 depth-two trees underfit
# (RMSE 0.626 versus 0.509 for the forest) and produced a flat focal interval,
# whereas 1,000 depth-three trees reached comparable RMSE (0.500) and resolved
# Longitude. None of the confirmatory seeds was inspected for this choice.

source(file.path("exploration", "s35_california_longitude_robustness.R"))

confirmatory_success_criteria <- function(repetitions) {
  list(
    minimum_feasible_splits = ceiling(0.95 * repetitions),
    minimum_sign_flip_fraction = 0.80,
    minimum_median_ess = 100,
    minimum_median_scaled_distance = 0.05
  )
}

fit_confirmatory_models <- function(training, response, seed) {
  set.seed(seed + 2L)
  forest <- randomForest::randomForest(
    x = training,
    y = response,
    ntree = 100L,
    mtry = max(1L, floor(sqrt(ncol(training)))),
    nodesize = 5L
  )
  set.seed(seed + 3L)
  boosting <- gbm::gbm.fit(
    x = training,
    y = response,
    distribution = "gaussian",
    n.trees = 1000L,
    interaction.depth = 3L,
    n.minobsinnode = 10L,
    shrinkage = 0.05,
    bag.fraction = 0.8,
    nTrain = nrow(training),
    keep.data = FALSE,
    verbose = FALSE
  )
  list(
    "Random forest" = function(newdata) {
      as.numeric(stats::predict(forest, newdata = newdata))
    },
    "Gradient boosting" = function(newdata) {
      as.numeric(stats::predict(
        boosting, newdata = newdata, n.trees = 1000L,
        type = "response"
      ))
    }
  )
}

empirical_reference_covariance <- function(predictions, weights) {
  predictions <- as.matrix(predictions)
  weights <- as.numeric(weights)
  if (nrow(predictions) != length(weights) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("weights must be finite, nonnegative, and match prediction rows.",
         call. = FALSE)
  }
  weights <- weights / sum(weights)
  density_ratio <- nrow(predictions) * weights
  ordinary_pd <- colMeans(predictions)
  centered_predictions <- sweep(predictions, 2L, ordinary_pd, "-")
  covariance_path <- colMeans((density_ratio - 1) * centered_predictions)
  direct_difference <- as.numeric(crossprod(weights, predictions)) - ordinary_pd
  list(
    covariance_path = covariance_path,
    direct_difference = direct_difference,
    identity_error = abs(covariance_path - direct_difference)
  )
}

fit_confirmatory_split <- function(data, split_id, seed,
                                   epsilon_values = c(0.10, 0.15, 0.20),
                                   alpha = 0.10, k = 21L) {
  split <- pilot_split(nrow(data$predictors), seed)
  split$evaluation <- pilot_cap_evaluation(
    split$evaluation, maximum = 5000L, seed = seed + 1L
  )
  training <- data$predictors[split$training, , drop = FALSE]
  calibration <- data$predictors[split$calibration, , drop = FALSE]
  evaluation <- data$predictors[split$evaluation, , drop = FALSE]
  response_training <- data$response[split$training]
  response_evaluation <- data$response[split$evaluation]

  limits <- stats::quantile(
    training$Longitude, c(0.40, 0.60), names = FALSE, type = 8
  )
  grid <- seq(limits[[1]], limits[[2]], length.out = k)
  learners <- list(
    "Locally scaled residual" = pilot_validity_rule(
      training, calibration, "Longitude", alpha
    ),
    "Conditional Gaussian density" = pilot_conditional_density_rule(
      training, calibration, "Longitude", alpha
    )
  )
  validity_matrices <- lapply(learners, function(learner) {
    validity_matrix(evaluation, "Longitude", grid, learner)
  })
  reference_epsilon <- 0.15
  reference_fits <- lapply(validity_matrices, function(validity) {
    fit_relaxed_reference_precise(validity, epsilon = reference_epsilon)
  })
  context_names <- setdiff(names(evaluation), "Longitude")
  balance_rows <- lapply(names(reference_fits), function(learner) {
    weights <- reference_fits[[learner]]$weights
    original_mean <- vapply(
      evaluation[, context_names, drop = FALSE], mean, numeric(1)
    )
    original_sd <- vapply(
      evaluation[, context_names, drop = FALSE], stats::sd, numeric(1)
    )
    weighted_mean <- vapply(context_names, function(variable) {
      sum(weights * evaluation[[variable]])
    }, numeric(1))
    data.frame(
      split_id = split_id,
      seed = seed,
      learner = learner,
      epsilon = reference_epsilon,
      feature = context_names,
      original_mean = original_mean,
      weighted_mean = weighted_mean,
      original_sd = original_sd,
      standardized_mean_shift =
        (weighted_mean - original_mean) / original_sd,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  models <- fit_confirmatory_models(
    training, response_training, seed
  )
  model_analyses <- lapply(names(models), function(model_name) {
    predict_function <- models[[model_name]]
    observed_predictions <- predict_function(evaluation)
    prediction_iqr <- diff(stats::quantile(
      observed_predictions, c(0.25, 0.75), names = FALSE
    ))
    rmse <- sqrt(mean((response_evaluation - observed_predictions)^2))
    predictions <- prediction_matrix(
      evaluation, "Longitude", grid, predict_function
    )
    ale_breaks <- unique(stats::quantile(
      evaluation$Longitude, probs = seq(0, 1, length.out = 21L),
      names = FALSE, type = 8
    ))
    ale_estimate <- estimate_ale(
      evaluation, "Longitude", ale_breaks, predict_function
    )
    ale_on_grid <- stats::approx(
      ale_estimate$z, ale_estimate$ale, xout = grid, rule = 2
    )$y

    learner_analyses <- lapply(names(validity_matrices), function(learner) {
      validity <- validity_matrices[[learner]]
      analysis <- analyze_validity_matrix(
        validity, predictions, prediction_iqr, epsilon_values,
        split_id, seed, learner, rmse, ale = ale_on_grid
      )
      analysis$summary$model <- model_name
      covariance <- empirical_reference_covariance(
        predictions, reference_fits[[learner]]$weights
      )
      covariance_rows <- data.frame(
        split_id = split_id,
        seed = seed,
        model = model_name,
        learner = learner,
        epsilon = reference_epsilon,
        grid_fraction = seq(0, 1, length.out = k),
        z = grid,
        covariance = covariance$covariance_path,
        centered_covariance = covariance$covariance_path -
          mean(covariance$covariance_path),
        direct_difference = covariance$direct_difference,
        centered_direct_difference = covariance$direct_difference -
          mean(covariance$direct_difference),
        identity_error = covariance$identity_error,
        stringsAsFactors = FALSE
      )
      curve_rows <- lapply(seq_along(epsilon_values), function(index) {
        epsilon <- epsilon_values[[index]]
        data.frame(
          split_id = split_id,
          seed = seed,
          model = model_name,
          learner = learner,
          epsilon = epsilon,
          grid_fraction = seq(0, 1, length.out = k),
          z = grid,
          ordinary_pd = colMeans(predictions),
          ale = analysis$ale,
          rcpd = analysis$rcpd[[index]],
          hard_cspd = analysis$hard,
          stringsAsFactors = FALSE
        )
      })
      list(
        summary = analysis$summary,
        curves = do.call(rbind, curve_rows),
        covariance = covariance_rows
      )
    })
    list(
      summary = do.call(rbind, lapply(learner_analyses, `[[`, "summary")),
      curves = do.call(rbind, lapply(learner_analyses, `[[`, "curves")),
      covariance = do.call(
        rbind, lapply(learner_analyses, `[[`, "covariance")
      )
    )
  })
  list(
    summary = do.call(rbind, lapply(model_analyses, `[[`, "summary")),
    curves = do.call(rbind, lapply(model_analyses, `[[`, "curves")),
    balance = do.call(rbind, balance_rows),
    covariance = do.call(rbind, lapply(model_analyses, `[[`, "covariance"))
  )
}

summarize_confirmatory_matrix <- function(results, repetitions) {
  criteria <- confirmatory_success_criteria(repetitions)
  groups <- split(results, interaction(
    results$model, results$learner, results$epsilon, drop = TRUE
  ))
  rows <- lapply(groups, function(group) {
    feasible <- group$rcpd_feasible &
      is.finite(group$rcpd_endpoint_contrast)
    usable <- group[feasible, , drop = FALSE]
    data.frame(
      model = group$model[[1]],
      learner = group$learner[[1]],
      epsilon = group$epsilon[[1]],
      feasible_splits = sum(feasible),
      sign_flip_fraction = mean(usable$pd_rcpd_sign_flip),
      rcpd_hard_sign_agreement_fraction =
        mean(usable$rcpd_hard_sign_agreement),
      median_pd_endpoint = stats::median(usable$pd_endpoint_contrast),
      median_ale_endpoint = stats::median(usable$ale_endpoint_contrast),
      median_rcpd_endpoint = stats::median(usable$rcpd_endpoint_contrast),
      median_hard_endpoint = stats::median(usable$hard_endpoint_contrast),
      median_pd_endpoint_normalized = stats::median(
        usable$pd_endpoint_over_prediction_iqr
      ),
      median_ale_endpoint_normalized = stats::median(
        usable$ale_endpoint_over_prediction_iqr
      ),
      median_rcpd_endpoint_normalized = stats::median(
        usable$rcpd_endpoint_over_prediction_iqr
      ),
      median_hard_endpoint_normalized = stats::median(
        usable$hard_endpoint_over_prediction_iqr
      ),
      median_gamma = stats::median(usable$gamma),
      median_ess = stats::median(usable$rcpd_ess),
      q10_ess = unname(stats::quantile(usable$rcpd_ess, 0.10, type = 8)),
      median_ess_fraction = stats::median(usable$rcpd_ess_fraction),
      median_scaled_distance = stats::median(
        usable$rcpd_distance_over_prediction_iqr
      ),
      median_rmse = stats::median(usable$rmse),
      max_constraint_violation = max(usable$rcpd_max_violation),
      max_kkt_residual = max(usable$rcpd_kkt_residual),
      passes_feasibility = sum(feasible) >=
        criteria$minimum_feasible_splits,
      passes_sign_flip = mean(usable$pd_rcpd_sign_flip) >=
        criteria$minimum_sign_flip_fraction,
      passes_ess = stats::median(usable$rcpd_ess) >=
        criteria$minimum_median_ess,
      passes_distance = stats::median(
        usable$rcpd_distance_over_prediction_iqr
      ) >= criteria$minimum_median_scaled_distance,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary$overall_pass <- with(
    summary,
    passes_feasibility & passes_sign_flip & passes_ess & passes_distance
  )
  summary$shape_sensitivity_pass <- with(
    summary,
    passes_feasibility & passes_ess & passes_distance
  )
  summary[order(summary$learner, summary$epsilon, summary$model), , drop = FALSE]
}

summarize_model_wide_gate <- function(summary) {
  settings <- split(summary, interaction(
    summary$learner, summary$epsilon, drop = TRUE
  ))
  rows <- lapply(settings, function(setting) {
    data.frame(
      learner = setting$learner[[1]],
      epsilon = setting$epsilon[[1]],
      model_classes = nrow(setting),
      passing_model_classes = sum(setting$overall_pass),
      all_models_pass = all(setting$overall_pass),
      shape_sensitive_model_classes = sum(setting$shape_sensitivity_pass),
      all_models_shape_sensitive = all(setting$shape_sensitivity_pass),
      minimum_sign_flip_fraction = min(setting$sign_flip_fraction),
      minimum_median_ess = min(setting$median_ess),
      minimum_median_scaled_distance = min(setting$median_scaled_distance),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result[order(result$learner, result$epsilon), , drop = FALSE]
}

plot_confirmatory_matrix <- function(results, path) {
  pd_rows <- unique(results[, c(
    "split_id", "seed", "model", "learner", "pd_endpoint_contrast"
  )])
  names(pd_rows)[names(pd_rows) == "pd_endpoint_contrast"] <- "contrast"
  pd_rows$setting <- "Ordinary PD"
  rcpd_rows <- results[, c(
    "split_id", "seed", "model", "learner", "epsilon",
    "rcpd_endpoint_contrast"
  )]
  names(rcpd_rows)[names(rcpd_rows) == "rcpd_endpoint_contrast"] <- "contrast"
  rcpd_rows$setting <- sprintf("RCPD, epsilon = %.2f", rcpd_rows$epsilon)
  plot_data <- rbind(
    pd_rows[, c(
      "split_id", "seed", "model", "learner", "setting", "contrast"
    )],
    rcpd_rows[, c(
      "split_id", "seed", "model", "learner", "setting", "contrast"
    )]
  )
  plot_data$setting <- factor(
    plot_data$setting,
    levels = c(
      "Ordinary PD", "RCPD, epsilon = 0.10",
      "RCPD, epsilon = 0.15", "RCPD, epsilon = 0.20"
    )
  )
  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(setting, contrast, color = setting)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, color = pdp_palette[["mid_grey"]], linetype = "dotted"
    ) +
    ggplot2::geom_boxplot(
      width = 0.55, outlier.shape = NA, color = pdp_palette[["mid_grey"]],
      fill = NA, linewidth = 0.45
    ) +
    ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.15, alpha = 0.65) +
    ggplot2::facet_grid(model ~ learner) +
    ggplot2::scale_color_manual(
      values = c(
        "Ordinary PD" = pdp_palette[["ink"]],
        "RCPD, epsilon = 0.10" = pdp_palette[["blue"]],
        "RCPD, epsilon = 0.15" = pdp_palette[["green"]],
        "RCPD, epsilon = 0.20" = pdp_palette[["orange"]]
      ),
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = "Longitude endpoint contrast") +
    theme_pdp(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 22, hjust = 1)
    )
  save_pdp_plot(figure, path, width = 10.0, height = 7.0)
}

run_longitude_confirmatory <- function(
    repetitions = 30L,
    base_seed = 20261101L,
    epsilon_values = c(0.10, 0.15, 0.20)) {
  data <- pilot_california_data()
  seeds <- base_seed + seq_len(repetitions) - 1L
  analyses <- lapply(seq_along(seeds), function(i) {
    message(sprintf("Split %d/%d (seed %d)", i, repetitions, seeds[[i]]))
    fit_confirmatory_split(
      data, i, seeds[[i]], epsilon_values = epsilon_values
    )
  })
  results <- do.call(rbind, lapply(analyses, `[[`, "summary"))
  curves <- do.call(rbind, lapply(analyses, `[[`, "curves"))
  summary <- summarize_confirmatory_matrix(results, repetitions)
  gate <- summarize_model_wide_gate(summary)

  output_directory <- file.path("results", "pilot")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    results,
    file.path(output_directory, "california-longitude-confirmatory.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(
      output_directory, "california-longitude-confirmatory-summary.csv"
    ),
    row.names = FALSE
  )
  utils::write.csv(
    gate,
    file.path(output_directory, "california-longitude-confirmatory-gate.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    curves,
    file.path(output_directory, "california-longitude-confirmatory-curves.csv"),
    row.names = FALSE
  )
  saveRDS(
    list(
      results = results,
      summary = summary,
      gate = gate,
      curves = curves,
      design = list(
        repetitions = repetitions,
        seeds = seeds,
        model_classes = c("Random forest", "Gradient boosting"),
        learners = c(
          "Locally scaled residual", "Conditional Gaussian density"
        ),
        epsilon_values = epsilon_values,
        alpha = 0.10,
          grid_quantiles = c(0.40, 0.60),
          k = 21L,
          random_forest_trees = 100L,
        boosting_trees = 1000L,
        boosting_depth = 3L,
        boosting_shrinkage = 0.05,
        criteria = confirmatory_success_criteria(repetitions)
      )
    ),
    file.path(output_directory, "california-longitude-confirmatory.rds")
  )
  plot_confirmatory_matrix(
    results,
    file.path(output_directory, "california-longitude-confirmatory.png")
  )
  list(summary = summary, gate = gate)
}

if (sys.nframe() == 0L) {
  confirmatory <- run_longitude_confirmatory()
  cat("\nModel-specific results:\n")
  print(confirmatory$summary, digits = 4, row.names = FALSE)
  cat("\nModel-wide gate:\n")
  print(confirmatory$gate, digits = 4, row.names = FALSE)
}
