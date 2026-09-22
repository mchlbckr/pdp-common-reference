# Confirmatory California Housing application.
#
# The implementation is loaded in a private environment because the
# exploratory development scripts contain generic helper names that should not
# alter the remaining targets pipeline. The confirmatory design itself is fixed
# in s36_california_longitude_confirmatory.R.

s33_environment <- new.env(parent = globalenv())
for (script in c(
  "s33_real_data_screen.R",
  "s34_california_longitude_stability.R",
  "s35_california_longitude_robustness.R",
  "s36_california_longitude_confirmatory.R"
)) {
  sys.source(file.path("exploration", script), envir = s33_environment)
}

run_s33 <- function(repetitions = 30L, base_seed = 20261101L,
                    epsilon_values = c(0.10, 0.15, 0.20)) {
  data <- s33_environment$pilot_california_data()
  seeds <- base_seed + seq_len(repetitions) - 1L
  analyses <- lapply(seq_along(seeds), function(i) {
    message(sprintf(
      "California confirmatory split %d/%d (seed %d)",
      i, repetitions, seeds[[i]]
    ))
    s33_environment$fit_confirmatory_split(
      data, i, seeds[[i]], epsilon_values = epsilon_values
    )
  })
  results <- do.call(rbind, lapply(analyses, `[[`, "summary"))
  curves <- do.call(rbind, lapply(analyses, `[[`, "curves"))
  balance <- do.call(rbind, lapply(analyses, `[[`, "balance"))
  covariance <- do.call(rbind, lapply(analyses, `[[`, "covariance"))
  summary <- s33_environment$summarize_confirmatory_matrix(
    results, repetitions
  )
  gate <- s33_environment$summarize_model_wide_gate(summary)
  list(
    results = results,
    curves = curves,
    balance = balance,
    covariance = covariance,
    summary = summary,
    gate = gate,
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
      criteria = s33_environment$confirmatory_success_criteria(repetitions)
    )
  )
}

summarize_s33 <- function(simulation) simulation$summary

gate_s33 <- function(simulation) simulation$gate

s33_quantile_summary <- function(x) {
  values <- stats::quantile(
    x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE,
    names = FALSE, type = 8
  )
  c(q25 = values[[1]], median = values[[2]], q75 = values[[3]])
}

audit_s33 <- function(simulation) {
  # Support and projection diagnostics do not depend on the prediction model.
  # Select one model to avoid counting every split twice.
  results <- simulation$results
  results <- results[
    results$model == simulation$design$model_classes[[1]], , drop = FALSE
  ]
  measures <- c(
    "min_isc", "mean_isc", "gamma", "common_count", "epsilon_min",
    "observed_patterns", "median_pattern_size", "min_pattern_size",
    "singleton_pattern_fraction", "rcpd_min_coverage",
    "rcpd_mean_coverage", "rcpd_simultaneous_valid_mass", "rcpd_kl",
    "rcpd_max_weight", "rcpd_ess", "rcpd_ess_fraction",
    "rcpd_active_constraints", "rcpd_active_affine_dimension",
    "rcpd_active_licq"
  )
  groups <- split(results, interaction(
    results$learner, results$epsilon, drop = TRUE
  ))
  output <- do.call(rbind, lapply(groups, function(group) {
    row <- data.frame(
      learner = group$learner[[1]], epsilon = group$epsilon[[1]],
      stringsAsFactors = FALSE
    )
    for (measure in measures) {
      quantiles <- s33_quantile_summary(group[[measure]])
      row[[paste0(measure, "_q25")]] <- quantiles[["q25"]]
      row[[paste0(measure, "_median")]] <- quantiles[["median"]]
      row[[paste0(measure, "_q75")]] <- quantiles[["q75"]]
    }
    row$rcpd_active_licq_fraction <- mean(group$rcpd_active_licq, na.rm = TRUE)
    row
  }))
  output[order(output$learner, output$epsilon), , drop = FALSE]
}

reference_diagnostics_s33 <- function(simulation) {
  balance_groups <- split(
    simulation$balance,
    interaction(simulation$balance$learner, simulation$balance$feature,
                drop = TRUE)
  )
  balance <- do.call(rbind, lapply(balance_groups, function(group) {
    quantiles <- s33_quantile_summary(group$standardized_mean_shift)
    data.frame(
      learner = group$learner[[1]],
      feature = group$feature[[1]],
      q25 = quantiles[["q25"]],
      median = quantiles[["median"]],
      q75 = quantiles[["q75"]],
      stringsAsFactors = FALSE
    )
  }))
  covariance_groups <- split(
    simulation$covariance,
    interaction(
      simulation$covariance$model, simulation$covariance$learner,
      simulation$covariance$grid_fraction, drop = TRUE
    )
  )
  covariance <- do.call(rbind, lapply(covariance_groups, function(group) {
    quantiles <- s33_quantile_summary(group$centered_covariance)
    data.frame(
      model = group$model[[1]],
      learner = group$learner[[1]],
      grid_fraction = group$grid_fraction[[1]],
      z = stats::median(group$z),
      q25 = quantiles[["q25"]],
      median = quantiles[["median"]],
      q75 = quantiles[["q75"]],
      stringsAsFactors = FALSE
    )
  }))
  list(
    balance = balance[order(balance$feature, balance$learner), , drop = FALSE],
    covariance = covariance[order(
      covariance$model, covariance$learner, covariance$grid_fraction
    ), , drop = FALSE],
    max_identity_error = max(simulation$covariance$identity_error),
    epsilon = unique(simulation$covariance$epsilon)
  )
}

plot_s33 <- function(simulation, directory = "results") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  results <- simulation$results
  pd_rows <- unique(results[, c(
    "split_id", "seed", "model", "learner", "pd_endpoint_contrast"
  )])
  names(pd_rows)[names(pd_rows) == "pd_endpoint_contrast"] <- "contrast"
  pd_rows$setting <- "Ordinary PD"
  ale_rows <- unique(results[, c(
    "split_id", "seed", "model", "learner", "ale_endpoint_contrast"
  )])
  names(ale_rows)[names(ale_rows) == "ale_endpoint_contrast"] <- "contrast"
  ale_rows$setting <- "ALE"
  hard_rows <- unique(results[, c(
    "split_id", "seed", "model", "learner", "hard_endpoint_contrast"
  )])
  names(hard_rows)[names(hard_rows) == "hard_endpoint_contrast"] <- "contrast"
  hard_rows$setting <- "Hard CSPD"
  rcpd_rows <- results[, c(
    "split_id", "seed", "model", "learner", "epsilon",
    "rcpd_endpoint_contrast"
  )]
  names(rcpd_rows)[names(rcpd_rows) == "rcpd_endpoint_contrast"] <- "contrast"
  rcpd_rows$setting <- sprintf("RCPD %.2f", rcpd_rows$epsilon)
  keep <- c("split_id", "seed", "model", "learner", "setting", "contrast")
  plot_data <- rbind(
    pd_rows[, keep], ale_rows[, keep], hard_rows[, keep], rcpd_rows[, keep]
  )
  setting_levels <- c(
    "Ordinary PD", "ALE", "Hard CSPD", "RCPD 0.10", "RCPD 0.15",
    "RCPD 0.20"
  )
  plot_data$setting <- factor(plot_data$setting, levels = setting_levels)
  colors <- c(
    "Ordinary PD" = pdp_palette[["ink"]],
    "ALE" = pdp_palette[["gold"]],
    "Hard CSPD" = pdp_palette[["purple"]],
    "RCPD 0.10" = pdp_palette[["blue"]],
    "RCPD 0.15" = pdp_palette[["green"]],
    "RCPD 0.20" = pdp_palette[["orange"]]
  )
  combinations <- data.frame(
    model = rep(c("Random forest", "Gradient boosting"), each = 2L),
    learner = rep(
      c("Conditional Gaussian density", "Locally scaled residual"),
      times = 2L
    ),
    stringsAsFactors = FALSE
  )
  paths <- file.path(
    directory,
    sprintf("s33-california-longitude-application-%s.png", letters[1:4])
  )
  y_limits <- range(plot_data$contrast, na.rm = TRUE)
  for (index in seq_len(nrow(combinations))) {
    selected <- plot_data$model == combinations$model[[index]] &
      plot_data$learner == combinations$learner[[index]]
    panel_data <- plot_data[selected, , drop = FALSE]
    panel <- ggplot2::ggplot(
      panel_data, ggplot2::aes(setting, contrast, color = setting)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, color = pdp_palette[["mid_grey"]],
        linetype = "dotted"
      ) +
      ggplot2::geom_boxplot(
        width = 0.55, outlier.shape = NA,
        color = pdp_palette[["mid_grey"]], fill = NA, linewidth = 0.45
      ) +
      ggplot2::geom_jitter(
        width = 0.10, height = 0, size = 1.05, alpha = 0.62
      ) +
      ggplot2::scale_color_manual(values = colors, drop = FALSE) +
      ggplot2::scale_x_discrete(
        labels = c("PD", "ALE", "Hard", "RCPD\n.10", "RCPD\n.15", "RCPD\n.20")
      ) +
      ggplot2::coord_cartesian(ylim = y_limits) +
      ggplot2::labs(x = NULL, y = "Longitude endpoint contrast") +
      theme_pdp(base_size = 9) +
      ggplot2::theme(legend.position = "none")
    save_pdp_plot(panel, paths[[index]], width = 5.0, height = 3.6)
  }
  paths
}

plot_s33_reference_diagnostics <- function(simulation, directory = "results") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- reference_diagnostics_s33(simulation)
  balance <- diagnostics$balance
  covariance <- diagnostics$covariance
  learner_labels <- c(
    "Conditional Gaussian density" = "Density",
    "Locally scaled residual" = "Scaled residual"
  )
  model_colors <- c(
    "Random forest" = pdp_palette[["blue"]],
    "Gradient boosting" = pdp_palette[["orange"]]
  )
  learner_colors <- c(
    "Conditional Gaussian density" = pdp_palette[["blue"]],
    "Locally scaled residual" = pdp_palette[["orange"]]
  )
  feature_labels <- c(
    AveBedrms = "Average bedrooms", AveOccup = "Average occupancy",
    AveRooms = "Average rooms", HouseAge = "House age",
    Latitude = "Latitude", MedInc = "Median income",
    Population = "Population"
  )
  feature_order <- names(sort(tapply(
    abs(balance$median), balance$feature, max
  )))
  balance$feature_index <- match(balance$feature, feature_order)
  balance$vertical_offset <- ifelse(
    balance$learner == "Conditional Gaussian density", -0.14, 0.14
  )
  balance$y <- balance$feature_index + balance$vertical_offset
  panel_a <- ggplot2::ggplot(balance) +
    ggplot2::geom_vline(
      xintercept = 0, color = pdp_palette[["mid_grey"]],
      linetype = "dotted"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = q25, xend = q75, y = y, yend = y, color = learner
      ), linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = median, y = y, color = learner, shape = learner),
      size = 2.1
    ) +
    ggplot2::scale_color_manual(
      values = learner_colors, labels = learner_labels
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Conditional Gaussian density" = 16,
        "Locally scaled residual" = 17
      ), labels = learner_labels
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(feature_order), labels = feature_labels[feature_order],
      expand = ggplot2::expansion(add = 0.45)
    ) +
    ggplot2::labs(
      x = "Standardized mean shift under the projected reference",
      y = NULL, color = "Validity rule", shape = "Validity rule"
    ) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  covariance$group <- interaction(
    covariance$model, covariance$learner, drop = TRUE
  )
  panel_b <- ggplot2::ggplot(
    covariance,
    ggplot2::aes(
      x = z, y = median, color = model, linetype = learner,
      group = group
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, color = pdp_palette[["mid_grey"]],
      linetype = "dotted"
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q25, ymax = q75, fill = model),
      color = NA, alpha = 0.10, show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(
      values = model_colors,
      labels = c("Gradient boosting" = "Boosting", "Random forest" = "Forest")
    ) +
    ggplot2::scale_fill_manual(values = model_colors) +
    ggplot2::scale_linetype_manual(
      values = c(
        "Conditional Gaussian density" = "solid",
        "Locally scaled residual" = "longdash"
      ), labels = learner_labels
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, nrow = 1),
      linetype = ggplot2::guide_legend(order = 2, nrow = 1)
    ) +
    ggplot2::labs(
      x = "Longitude", y = "Centered covariance path",
      color = "Prediction model", linetype = "Validity rule"
    ) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.box.just = "left"
    )

  paths <- file.path(
    directory,
    c(
      "s33-california-reference-diagnostics-a.png",
      "s33-california-reference-diagnostics-b.png"
    )
  )
  save_pdp_plot(panel_a, paths[[1]], width = 5.2, height = 3.9)
  save_pdp_plot(panel_b, paths[[2]], width = 5.2, height = 3.9)
  paths
}
