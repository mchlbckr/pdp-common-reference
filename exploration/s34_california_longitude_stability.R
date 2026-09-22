# Split-stability check for the California-housing Longitude pilot.
#
# This is a deliberately separate, exploratory analysis. It does not feed the
# targets pipeline, paper, or supplement. The design is fixed before execution:
# 20 train/calibration/evaluation splits, the central 20% of training Longitude,
# K = 21, alpha = epsilon = 0.10, and the same random-forest specification used
# in s33_real_data_screen.R.

source(file.path("exploration", "s33_real_data_screen.R"))

stability_success_criteria <- function() {
  list(
    minimum_feasible_splits = 19L,
    minimum_sign_flip_fraction = 0.80,
    minimum_median_ess = 100,
    minimum_median_scaled_distance = 0.05
  )
}

# The production helper uses the default L-BFGS-B stopping rule and then checks
# feasibility at a tighter scale. In this repeated-split diagnostic that rule
# occasionally stops 2e-5--6e-5 before the boundary. The locally defined
# refinement keeps the estimand unchanged and uses a stricter objective stop.
fit_relaxed_reference_precise <- function(validity, epsilon = 0.10,
                                          tolerance = 1e-6) {
  h <- unclass(validity) * 1
  target <- 1 - epsilon
  objective <- function(lambda) {
    eta <- as.numeric(h %*% lambda)
    log_mean_exp(eta) - target * sum(lambda)
  }
  gradient <- function(lambda) {
    eta <- as.numeric(h %*% lambda)
    weights <- exp(eta - max(eta))
    weights <- weights / sum(weights)
    as.numeric(crossprod(weights, h)) - target
  }
  fit <- stats::optim(
    par = rep(0, ncol(h)), fn = objective, gr = gradient,
    method = "L-BFGS-B", lower = rep(0, ncol(h)),
    upper = rep(50, ncol(h)),
    control = list(maxit = 20000L, factr = 10, pgtol = 1e-12)
  )
  eta <- as.numeric(h %*% fit$par)
  weights <- exp(eta - max(eta))
  weights <- weights / sum(weights)
  coverage <- as.numeric(crossprod(weights, h))
  max_violation <- max(0, target - min(coverage))
  gradient_at_solution <- coverage - target
  active <- fit$par > sqrt(.Machine$double.eps)
  kkt_residual <- max(c(
    abs(gradient_at_solution[active]),
    pmax(0, -gradient_at_solution[!active])
  ))
  validity_diagnostics <- reference_validity_diagnostics(validity, weights)
  pattern_diagnostics <- validity_pattern_diagnostics(
    validity, weights = weights, multipliers = fit$par,
    tolerance = tolerance
  )
  list(
    weights = weights,
    coverage = coverage,
    multipliers = fit$par,
    convergence = fit$convergence,
    max_violation = max_violation,
    kkt_residual = kkt_residual,
    feasible = max_violation <= tolerance && kkt_residual <= tolerance,
    effective_sample_size = 1 / sum(weights^2),
    kl = sum(weights * log(weights * nrow(h))),
    min_coverage = validity_diagnostics$min_coverage,
    mean_coverage = validity_diagnostics$mean_coverage,
    simultaneous_valid_mass = validity_diagnostics$simultaneous_valid_mass,
    max_weight = max(weights),
    epsilon_min = pattern_diagnostics$epsilon_min,
    observed_patterns = pattern_diagnostics$observed_patterns,
    median_pattern_size = pattern_diagnostics$median_pattern_size,
    min_pattern_size = pattern_diagnostics$min_pattern_size,
    singleton_pattern_fraction = pattern_diagnostics$singleton_pattern_fraction,
    active_constraints = pattern_diagnostics$active_constraints,
    active_affine_dimension = pattern_diagnostics$active_affine_dimension,
    active_licq = pattern_diagnostics$active_licq
  )
}

fit_longitude_split <- function(data, split_id, seed, alpha = 0.10,
                                epsilon = 0.10, k = 21L) {
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

  result <- pilot_analyze_feature(
    dataset = "California housing",
    focal = "Longitude",
    training = training,
    calibration = calibration,
    evaluation = evaluation,
    predict_function = predict_function,
    observed_prediction_iqr = prediction_iqr,
    alpha = alpha,
    epsilon = epsilon,
    grid_quantiles = c(0.40, 0.60),
    k = k
  )

  precise_fit <- fit_relaxed_reference_precise(
    result$validity, epsilon = epsilon
  )
  relaxed <- as.numeric(crossprod(precise_fit$weights, result$predictions))
  pdp <- colMeans(result$predictions)
  rcpd_rows <- result$curves$method == "RCPD"
  result$curves$value[rcpd_rows] <- relaxed
  result$curves$centered[rcpd_rows] <- relaxed - mean(relaxed)
  result$curves$scaled_centered[rcpd_rows] <-
    (relaxed - mean(relaxed)) / max(prediction_iqr, sqrt(.Machine$double.eps))

  summary <- result$summary
  summary$rcpd_feasible <- precise_fit$feasible
  summary$rcpd_ess <- precise_fit$effective_sample_size
  summary$rcpd_kl <- precise_fit$kl
  summary$rcpd_max_violation <- precise_fit$max_violation
  summary$rcpd_distance <- pilot_safe_distance(relaxed, pdp)
  summary$rcpd_distance_over_pd_amplitude <-
    summary$rcpd_distance / max(summary$pd_amplitude, sqrt(.Machine$double.eps))
  summary$rcpd_distance_over_prediction_iqr <-
    summary$rcpd_distance / max(prediction_iqr, sqrt(.Machine$double.eps))
  summary$rcpd_endpoint_contrast <- tail(relaxed, 1) - head(relaxed, 1)
  summary$split_id <- split_id
  summary$seed <- seed
  summary$n_training <- nrow(training)
  summary$n_calibration <- nrow(calibration)
  summary$n_evaluation <- nrow(evaluation)
  summary$rmse <- rmse
  summary$common_count <- sum(rowSums(result$validity) == ncol(result$validity))
  summary$rcpd_ess_fraction <- summary$rcpd_ess / nrow(evaluation)
  summary$pd_rcpd_sign_flip <- with(
    summary,
    is.finite(pd_endpoint_contrast) & is.finite(rcpd_endpoint_contrast) &
      sign(pd_endpoint_contrast) != sign(rcpd_endpoint_contrast)
  )
  summary$rcpd_hard_sign_agreement <- with(
    summary,
    is.finite(rcpd_endpoint_contrast) & is.finite(hard_endpoint_contrast) &
      sign(rcpd_endpoint_contrast) == sign(hard_endpoint_contrast)
  )

  curves <- result$curves
  curves$split_id <- split_id
  curves$seed <- seed
  curves$grid_fraction <- rep(seq(0, 1, length.out = k), times = 4L)

  list(summary = summary, curves = curves)
}

summarize_stability <- function(results) {
  criteria <- stability_success_criteria()
  feasible <- results$rcpd_feasible & is.finite(results$rcpd_endpoint_contrast)
  feasible_results <- results[feasible, , drop = FALSE]
  quantile_text <- function(x) {
    values <- stats::quantile(x, c(0.10, 0.50, 0.90), na.rm = TRUE,
                              names = FALSE, type = 8)
    stats::setNames(values, c("q10", "median", "q90"))
  }

  values <- list(
    repetitions = nrow(results),
    feasible_splits = sum(feasible),
    sign_flip_fraction = mean(feasible_results$pd_rcpd_sign_flip),
    rcpd_hard_sign_agreement_fraction =
      mean(feasible_results$rcpd_hard_sign_agreement),
    pd_endpoint_contrast = quantile_text(
      feasible_results$pd_endpoint_contrast
    ),
    rcpd_endpoint_contrast = quantile_text(
      feasible_results$rcpd_endpoint_contrast
    ),
    hard_endpoint_contrast = quantile_text(
      feasible_results$hard_endpoint_contrast
    ),
    gamma = quantile_text(feasible_results$gamma),
    rcpd_ess = quantile_text(feasible_results$rcpd_ess),
    rcpd_ess_fraction = quantile_text(
      feasible_results$rcpd_ess_fraction
    ),
    scaled_distance = quantile_text(
      feasible_results$rcpd_distance_over_prediction_iqr
    ),
    rmse = quantile_text(feasible_results$rmse)
  )
  values$criteria <- criteria
  values$passes <- c(
    feasible_splits = values$feasible_splits >=
      criteria$minimum_feasible_splits,
    sign_flip_fraction = values$sign_flip_fraction >=
      criteria$minimum_sign_flip_fraction,
    median_ess = values$rcpd_ess[["median"]] >=
      criteria$minimum_median_ess,
    median_scaled_distance = values$scaled_distance[["median"]] >=
      criteria$minimum_median_scaled_distance
  )
  values$overall_pass <- all(values$passes)
  values
}

plot_stability_results <- function(results, path) {
  contrasts <- rbind(
    data.frame(
      split_id = results$split_id,
      method = "Ordinary PD",
      contrast = results$pd_endpoint_contrast
    ),
    data.frame(
      split_id = results$split_id,
      method = "RCPD",
      contrast = results$rcpd_endpoint_contrast
    ),
    data.frame(
      split_id = results$split_id,
      method = "Hard CSPD",
      contrast = results$hard_endpoint_contrast
    )
  )
  contrasts$method <- factor(
    contrasts$method,
    levels = c("Ordinary PD", "RCPD", "Hard CSPD")
  )
  colors <- c(
    "Ordinary PD" = pdp_palette[["ink"]],
    "RCPD" = pdp_palette[["blue"]],
    "Hard CSPD" = pdp_palette[["orange"]]
  )
  figure <- ggplot2::ggplot(
    contrasts,
    ggplot2::aes(method, contrast, color = method)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, color = pdp_palette[["mid_grey"]], linetype = "dotted"
    ) +
    ggplot2::geom_boxplot(
      width = 0.50, outlier.shape = NA, color = pdp_palette[["mid_grey"]],
      fill = NA, linewidth = 0.45
    ) +
    ggplot2::geom_jitter(width = 0.10, height = 0, size = 1.8, alpha = 0.75) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "Endpoint contrast over central 20% of Longitude"
    ) +
    theme_pdp(base_size = 10) +
    ggplot2::theme(legend.position = "none")
  save_pdp_plot(figure, path, width = 6.5, height = 4.2)
}

run_longitude_stability <- function(repetitions = 20L,
                                    base_seed = 20261001L) {
  data <- pilot_california_data()
  seeds <- base_seed + seq_len(repetitions) - 1L
  analyses <- lapply(seq_along(seeds), function(i) {
    message(sprintf("Split %d/%d (seed %d)", i, repetitions, seeds[[i]]))
    fit_longitude_split(data, i, seeds[[i]])
  })
  results <- do.call(rbind, lapply(analyses, `[[`, "summary"))
  curves <- do.call(rbind, lapply(analyses, `[[`, "curves"))
  stability <- summarize_stability(results)

  output_directory <- file.path("results", "pilot")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    results,
    file.path(output_directory, "california-longitude-stability.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    curves,
    file.path(output_directory, "california-longitude-stability-curves.csv"),
    row.names = FALSE
  )
  saveRDS(
    list(
      results = results,
      curves = curves,
      stability = stability,
      design = list(
        repetitions = repetitions,
        seeds = seeds,
        grid_quantiles = c(0.40, 0.60),
        k = 21L,
        alpha = 0.10,
        epsilon = 0.10,
        random_forest_trees = 100L
      )
    ),
    file.path(output_directory, "california-longitude-stability.rds")
  )
  plot_stability_results(
    results,
    file.path(output_directory, "california-longitude-stability.png")
  )
  stability
}

if (sys.nframe() == 0L) {
  stability <- run_longitude_stability()
  print(stability)
}
