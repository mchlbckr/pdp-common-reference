# Robustness matrix for the California-housing Longitude pilot.
#
# The design is fixed before execution: the same 20 splits used in s34, two
# independently calibrated validity learners, epsilon in {0.10, 0.15, 0.20},
# the central 20% of training-set Longitude, K = 21, and the same 100-tree
# random forest. This exploratory script does not feed the paper pipeline.

source(file.path("exploration", "s34_california_longitude_stability.R"))

pilot_conditional_density_rule <- function(training, calibration, focal,
                                           alpha = 0.10) {
  standardizer <- pilot_standardizer(training)
  fit <- suppressWarnings(fit_conditional_density_validity(
    pilot_standardize(training, standardizer),
    pilot_standardize(calibration, standardizer),
    focal,
    alpha = alpha
  ))
  function(newdata) {
    suppressWarnings(predict_conditional_density_validity(
      fit,
      pilot_standardize(newdata, standardizer)
    ))
  }
}

analyze_validity_matrix <- function(validity, predictions, prediction_iqr,
                                    epsilon_values, split_id, seed, learner,
                                    rmse, ale = NULL) {
  pdp <- colMeans(predictions)
  if (is.null(ale)) ale <- rep(NA_real_, ncol(predictions))
  if (length(ale) != ncol(predictions)) {
    stop("ale must be NULL or match the prediction grid.", call. = FALSE)
  }
  common <- rowSums(validity) == ncol(validity)
  hard <- if (any(common)) {
    colMeans(predictions[common, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(predictions))
  }
  hard_endpoint <- if (all(is.finite(hard))) {
    tail(hard, 1) - head(hard, 1)
  } else {
    NA_real_
  }
  pd_endpoint <- tail(pdp, 1) - head(pdp, 1)
  ale_endpoint <- if (all(is.finite(ale))) {
    tail(ale, 1) - head(ale, 1)
  } else {
    NA_real_
  }
  pd_amplitude <- diff(range(pdp - mean(pdp)))

  rcpd_curves <- vector("list", length(epsilon_values))
  rows <- lapply(seq_along(epsilon_values), function(index) {
    epsilon <- epsilon_values[[index]]
    fit <- fit_relaxed_reference_precise(validity, epsilon = epsilon)
    rcpd <- as.numeric(crossprod(fit$weights, predictions))
    rcpd_curves[[index]] <<- rcpd
    rcpd_endpoint <- tail(rcpd, 1) - head(rcpd, 1)
    distance <- pilot_safe_distance(rcpd, pdp)
    data.frame(
      split_id = split_id,
      seed = seed,
      learner = learner,
      epsilon = epsilon,
      rcpd_feasible = fit$feasible,
      rcpd_optimizer_convergence = fit$convergence,
      rcpd_max_violation = fit$max_violation,
      rcpd_kkt_residual = fit$kkt_residual,
      min_isc = min(colMeans(validity)),
      mean_isc = mean(validity),
      gamma = mean(common),
      common_count = sum(common),
      epsilon_min = fit$epsilon_min,
      observed_patterns = fit$observed_patterns,
      median_pattern_size = fit$median_pattern_size,
      min_pattern_size = fit$min_pattern_size,
      singleton_pattern_fraction = fit$singleton_pattern_fraction,
      rcpd_min_coverage = fit$min_coverage,
      rcpd_mean_coverage = fit$mean_coverage,
      rcpd_simultaneous_valid_mass = fit$simultaneous_valid_mass,
      rcpd_max_weight = fit$max_weight,
      rcpd_active_constraints = fit$active_constraints,
      rcpd_active_affine_dimension = fit$active_affine_dimension,
      rcpd_active_licq = as.integer(fit$active_licq),
      rcpd_ess = fit$effective_sample_size,
      rcpd_ess_fraction = fit$effective_sample_size / nrow(validity),
      rcpd_kl = fit$kl,
      pd_endpoint_contrast = pd_endpoint,
      ale_endpoint_contrast = ale_endpoint,
      rcpd_endpoint_contrast = rcpd_endpoint,
      hard_endpoint_contrast = hard_endpoint,
      pd_endpoint_over_prediction_iqr = pd_endpoint /
        max(prediction_iqr, sqrt(.Machine$double.eps)),
      ale_endpoint_over_prediction_iqr = ale_endpoint /
        max(prediction_iqr, sqrt(.Machine$double.eps)),
      rcpd_endpoint_over_prediction_iqr = rcpd_endpoint /
        max(prediction_iqr, sqrt(.Machine$double.eps)),
      hard_endpoint_over_prediction_iqr = hard_endpoint /
        max(prediction_iqr, sqrt(.Machine$double.eps)),
      pd_rcpd_sign_flip = sign(pd_endpoint) != sign(rcpd_endpoint),
      rcpd_hard_sign_agreement = is.finite(hard_endpoint) &&
        sign(rcpd_endpoint) == sign(hard_endpoint),
      rcpd_distance = distance,
      rcpd_distance_over_prediction_iqr = distance /
        max(prediction_iqr, sqrt(.Machine$double.eps)),
      rcpd_distance_over_pd_amplitude = distance /
        max(pd_amplitude, sqrt(.Machine$double.eps)),
      rmse = rmse,
      stringsAsFactors = FALSE
    )
  })
  list(
    summary = do.call(rbind, rows),
    pdp = pdp,
    ale = ale,
    hard = hard,
    rcpd = rcpd_curves
  )
}

fit_robustness_split <- function(data, split_id, seed,
                                 epsilon_values = c(0.10, 0.15, 0.20),
                                 alpha = 0.10, k = 21L) {
  split <- pilot_split(nrow(data$predictors), seed)
  split$evaluation <- pilot_cap_evaluation(
    split$evaluation, maximum = 5000L, seed = seed + 1L
  )
  training <- data$predictors[split$training, , drop = FALSE]
  calibration <- data$predictors[split$calibration, , drop = FALSE]
  evaluation <- data$predictors[split$evaluation, , drop = FALSE]

  set.seed(seed + 2L)
  model <- randomForest::randomForest(
    x = training,
    y = data$response[split$training],
    ntree = 100L,
    mtry = max(1L, floor(sqrt(ncol(training)))),
    nodesize = 5L
  )
  predict_function <- function(newdata) {
    as.numeric(stats::predict(model, newdata = newdata))
  }
  observed_predictions <- predict_function(evaluation)
  prediction_iqr <- diff(stats::quantile(
    observed_predictions, c(0.25, 0.75), names = FALSE
  ))
  rmse <- sqrt(mean(
    (data$response[split$evaluation] - observed_predictions)^2
  ))
  limits <- stats::quantile(
    training$Longitude, c(0.40, 0.60), names = FALSE, type = 8
  )
  grid <- seq(limits[[1]], limits[[2]], length.out = k)
  predictions <- prediction_matrix(
    evaluation, "Longitude", grid, predict_function
  )

  learners <- list(
    "Locally scaled residual" = pilot_validity_rule(
      training, calibration, "Longitude", alpha
    ),
    "Conditional Gaussian density" = pilot_conditional_density_rule(
      training, calibration, "Longitude", alpha
    )
  )

  analyses <- lapply(names(learners), function(learner) {
    validity <- validity_matrix(
      evaluation, "Longitude", grid, learners[[learner]]
    )
    analysis <- analyze_validity_matrix(
      validity, predictions, prediction_iqr, epsilon_values,
      split_id, seed, learner, rmse
    )
    curve_rows <- lapply(seq_along(epsilon_values), function(index) {
      epsilon <- epsilon_values[[index]]
      data.frame(
        split_id = split_id,
        seed = seed,
        learner = learner,
        epsilon = epsilon,
        grid_fraction = seq(0, 1, length.out = k),
        z = grid,
        ordinary_pd = predictions |> colMeans(),
        rcpd = analysis$rcpd[[index]],
        hard_cspd = analysis$hard,
        stringsAsFactors = FALSE
      )
    })
    list(
      summary = analysis$summary,
      curves = do.call(rbind, curve_rows)
    )
  })
  list(
    summary = do.call(rbind, lapply(analyses, `[[`, "summary")),
    curves = do.call(rbind, lapply(analyses, `[[`, "curves"))
  )
}

summarize_robustness_matrix <- function(results) {
  criteria <- stability_success_criteria()
  groups <- split(results, interaction(
    results$learner, results$epsilon, drop = TRUE
  ))
  rows <- lapply(groups, function(group) {
    feasible <- group$rcpd_feasible &
      is.finite(group$rcpd_endpoint_contrast)
    usable <- group[feasible, , drop = FALSE]
    data.frame(
      learner = group$learner[[1]],
      epsilon = group$epsilon[[1]],
      feasible_splits = sum(feasible),
      sign_flip_fraction = mean(usable$pd_rcpd_sign_flip),
      rcpd_hard_sign_agreement_fraction =
        mean(usable$rcpd_hard_sign_agreement),
      median_pd_endpoint = stats::median(usable$pd_endpoint_contrast),
      median_rcpd_endpoint = stats::median(usable$rcpd_endpoint_contrast),
      median_hard_endpoint = stats::median(usable$hard_endpoint_contrast),
      median_gamma = stats::median(usable$gamma),
      median_ess = stats::median(usable$rcpd_ess),
      q10_ess = unname(stats::quantile(usable$rcpd_ess, 0.10, type = 8)),
      median_ess_fraction = stats::median(usable$rcpd_ess_fraction),
      median_scaled_distance = stats::median(
        usable$rcpd_distance_over_prediction_iqr
      ),
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
  summary[order(summary$learner, summary$epsilon), , drop = FALSE]
}

plot_robustness_matrix <- function(results, path) {
  pd_rows <- unique(results[, c(
    "split_id", "seed", "learner", "pd_endpoint_contrast"
  )])
  names(pd_rows)[names(pd_rows) == "pd_endpoint_contrast"] <- "contrast"
  pd_rows$setting <- "Ordinary PD"
  rcpd_rows <- results[, c(
    "split_id", "seed", "learner", "epsilon", "rcpd_endpoint_contrast"
  )]
  names(rcpd_rows)[names(rcpd_rows) == "rcpd_endpoint_contrast"] <- "contrast"
  rcpd_rows$setting <- sprintf("RCPD, epsilon = %.2f", rcpd_rows$epsilon)
  plot_data <- rbind(
    pd_rows[, c("split_id", "seed", "learner", "setting", "contrast")],
    rcpd_rows[, c("split_id", "seed", "learner", "setting", "contrast")]
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
    ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.4, alpha = 0.70) +
    ggplot2::facet_wrap(~ learner, ncol = 1) +
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
    theme_pdp(base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
  save_pdp_plot(figure, path, width = 8.0, height = 7.0)
}

run_longitude_robustness <- function(repetitions = 20L,
                                     base_seed = 20261001L,
                                     epsilon_values = c(0.10, 0.15, 0.20)) {
  data <- pilot_california_data()
  seeds <- base_seed + seq_len(repetitions) - 1L
  analyses <- lapply(seq_along(seeds), function(i) {
    message(sprintf("Split %d/%d (seed %d)", i, repetitions, seeds[[i]]))
    fit_robustness_split(
      data, i, seeds[[i]], epsilon_values = epsilon_values
    )
  })
  results <- do.call(rbind, lapply(analyses, `[[`, "summary"))
  curves <- do.call(rbind, lapply(analyses, `[[`, "curves"))
  summary <- summarize_robustness_matrix(results)

  output_directory <- file.path("results", "pilot")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    results,
    file.path(output_directory, "california-longitude-robustness.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(output_directory, "california-longitude-robustness-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    curves,
    file.path(output_directory, "california-longitude-robustness-curves.csv"),
    row.names = FALSE
  )
  saveRDS(
    list(
      results = results,
      summary = summary,
      curves = curves,
      design = list(
        repetitions = repetitions,
        seeds = seeds,
        learners = c(
          "Locally scaled residual", "Conditional Gaussian density"
        ),
        epsilon_values = epsilon_values,
        alpha = 0.10,
        grid_quantiles = c(0.40, 0.60),
        k = 21L,
        random_forest_trees = 100L
      )
    ),
    file.path(output_directory, "california-longitude-robustness.rds")
  )
  plot_robustness_matrix(
    results,
    file.path(output_directory, "california-longitude-robustness.png")
  )
  summary
}

if (sys.nframe() == 0L) {
  robustness <- run_longitude_robustness()
  print(robustness, digits = 4, row.names = FALSE)
}
