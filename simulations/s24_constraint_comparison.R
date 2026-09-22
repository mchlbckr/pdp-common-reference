# Simulation S24: rotating invalidity leaves no simultaneously valid context.
# Each context group is invalid at one grid value. The observed group mixture is
# deliberately imbalanced, so average validity can hide a severe local deficit.

s24_truth <- function(z) sin(2 * pi * (z - min(z)) / (max(z) - min(z)))

simulate_s24_data <- function(n, k = 8L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  probabilities <- exp(seq(-1.2, 1.2, length.out = k))
  probabilities <- probabilities / sum(probabilities)
  group <- sample.int(k, n, replace = TRUE, prob = probabilities)
  data.frame(x = rep(0, n), group = group)
}

s24_validity <- function(newdata) {
  newdata$group != as.integer(round(newdata$x))
}

s24_prediction <- function(k = 8L, beta = 1.2, off_support_strength = 3) {
  function(newdata) {
    z <- newdata$x
    group_scaled <- (newdata$group - mean(seq_len(k))) / k
    invalid <- newdata$group == as.integer(round(z))
    s24_truth(seq_len(k))[pmax(1L, pmin(k, as.integer(round(z))))] +
      beta * group_scaled + off_support_strength * invalid * group_scaled
  }
}

s24_one <- function(repetition, n_evaluation, k, epsilon, beta,
                    off_support_strength, seed) {
  data <- simulate_s24_data(n_evaluation, k, seed + repetition)
  grid <- seq_len(k)
  prediction <- s24_prediction(k, beta, off_support_strength)
  h <- validity_matrix(data, "x", grid, s24_validity)
  pointwise_reference <- fit_relaxed_reference(h, epsilon = epsilon,
                                               tolerance = 2e-5)
  average_reference <- fit_average_validity_reference(h, epsilon = epsilon)
  truth <- s24_truth(grid)

  pd <- estimate_pdp(data, "x", grid, prediction)$pdp
  trimmed <- estimate_pointwise_trimmed_pdp(
    data, "x", grid, s24_validity, prediction
  )$pointwise_trimmed_pdp
  pointwise <- estimate_weighted_pdp(
    data, "x", grid, prediction, pointwise_reference$weights
  )$weighted_pdp
  average <- estimate_weighted_pdp(
    data, "x", grid, prediction, average_reference$weights
  )$weighted_pdp
  curves <- list(
    PDP = pd, `Pointwise trim` = trimmed,
    `Pointwise constraints` = pointwise,
    `Average constraint` = average
  )
  truth_centered <- truth - mean(truth)
  truth_amplitude <- diff(range(truth_centered))
  metrics <- do.call(rbind, lapply(names(curves), function(method) {
    curve <- curves[[method]] - mean(curves[[method]])
    reference <- if (method == "Pointwise constraints") pointwise_reference else
      if (method == "Average constraint") average_reference else NULL
    data.frame(
      repetition = repetition, method = method,
      normalized_shape_error = max(abs(curve - truth_centered)) /
        truth_amplitude,
      min_coverage = if (is.null(reference)) NA_real_ else reference$min_coverage,
      mean_coverage = if (is.null(reference)) NA_real_ else reference$mean_coverage,
      simultaneous_valid_mass = if (is.null(reference)) 0 else
        reference$simultaneous_valid_mass,
      kl = if (is.null(reference)) NA_real_ else reference$kl,
      ess_fraction = if (is.null(reference)) NA_real_ else
        reference$effective_sample_size / n_evaluation,
      epsilon_min = pointwise_reference$epsilon_min,
      gamma = mean(rowSums(h) == ncol(h)),
      stringsAsFactors = FALSE
    )
  }))
  diagnostics <- data.frame(
    repetition = repetition,
    reference = c("Original P", "Pointwise constraints", "Average constraint",
                  "Intersection constraint"),
    feasible = c(TRUE, pointwise_reference$feasible,
                 average_reference$feasible, FALSE),
    min_coverage = c(min(colMeans(h)), pointwise_reference$min_coverage,
                     average_reference$min_coverage, NA_real_),
    mean_coverage = c(mean(h), pointwise_reference$mean_coverage,
                      average_reference$mean_coverage, NA_real_),
    simultaneous_valid_mass = c(0, pointwise_reference$simultaneous_valid_mass,
                                average_reference$simultaneous_valid_mass, NA_real_),
    kl = c(0, pointwise_reference$kl, average_reference$kl, NA_real_),
    ess_fraction = c(1, pointwise_reference$effective_sample_size / n_evaluation,
                     average_reference$effective_sample_size / n_evaluation,
                     NA_real_)
  )
  curve_data <- do.call(rbind, lapply(names(curves), function(method) {
    data.frame(z = grid, value = curves[[method]] - mean(curves[[method]]),
               method = method)
  }))
  list(metrics = metrics, diagnostics = diagnostics,
       curves = if (repetition == 1L) curve_data else NULL,
       truth = if (repetition == 1L) data.frame(z = grid, value = truth_centered) else NULL,
       observed_group_mass = if (repetition == 1L)
         data.frame(group = grid, mass = tabulate(data$group, nbins = k) / n_evaluation)
       else NULL)
}

run_s24 <- function(repetitions = 300L, n_evaluation = 1000L, k = 8L,
                    epsilon = 0.13, beta = 1.2, off_support_strength = 3,
                    seed = 20260926L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  if (epsilon < 1 / k) {
    stop("epsilon must be at least 1/k in the rotating-invalidity design.",
         call. = FALSE)
  }
  evaluate <- function(index) s24_one(
    index, n_evaluation, k, epsilon, beta, off_support_strength, seed
  )
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(repetitions), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(repetitions), evaluate)
  }
  list(
    metrics = do.call(rbind, lapply(pieces, `[[`, "metrics")),
    diagnostics = do.call(rbind, lapply(pieces, `[[`, "diagnostics")),
    curves = pieces[[1]]$curves, truth = pieces[[1]]$truth,
    observed_group_mass = pieces[[1]]$observed_group_mass,
    settings = list(repetitions = repetitions, n_evaluation = n_evaluation,
                    k = k, epsilon = epsilon, beta = beta,
                    off_support_strength = off_support_strength)
  )
}

summarize_s24 <- function(simulation) {
  metric_groups <- split(simulation$metrics, simulation$metrics$method)
  metrics <- do.call(rbind, lapply(metric_groups, function(part) {
    data.frame(
      method = part$method[[1]],
      normalized_shape_error = mean(part$normalized_shape_error),
      normalized_shape_error_mc_se = stats::sd(part$normalized_shape_error) /
        sqrt(nrow(part)),
      epsilon_min = mean(part$epsilon_min), gamma = mean(part$gamma),
      stringsAsFactors = FALSE
    )
  }))
  diagnostic_groups <- split(simulation$diagnostics,
                             simulation$diagnostics$reference)
  diagnostics <- do.call(rbind, lapply(diagnostic_groups, function(part) {
    data.frame(
      reference = part$reference[[1]], feasible_rate = mean(part$feasible),
      min_coverage = mean(part$min_coverage, na.rm = TRUE),
      mean_coverage = mean(part$mean_coverage, na.rm = TRUE),
      simultaneous_valid_mass = mean(part$simultaneous_valid_mass, na.rm = TRUE),
      kl = mean(part$kl, na.rm = TRUE),
      ess_fraction = mean(part$ess_fraction, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  list(metrics = metrics, diagnostics = diagnostics)
}

plot_s24 <- function(simulation, directory = "results") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  method_order <- c("Truth", "PDP", "Pointwise trim", "Pointwise constraints",
                    "Average constraint")
  colors <- c(
    Truth = pdp_palette[["ink"]], PDP = pdp_palette[["orange"]],
    `Pointwise trim` = pdp_palette[["purple"]],
    `Pointwise constraints` = pdp_palette[["blue"]],
    `Average constraint` = pdp_palette[["gold"]]
  )
  curves <- rbind(
    transform(simulation$truth, method = "Truth"), simulation$curves
  )
  curves$method <- factor(curves$method, levels = method_order)
  p_curves <- ggplot2::ggplot(
    curves, ggplot2::aes(z, value, color = method, linetype = method)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::scale_color_manual(
      values = colors,
      labels = c("Truth", "PDP", "Pointwise trim", "RCPD",
                 "Average-only projection")
    ) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash",
                                              "longdash", "dotted"),
                                    labels = c(
                                      "Truth", "PDP", "Pointwise trim",
                                      "RCPD", "Average-only projection"
                                    )) +
    ggplot2::scale_x_continuous(breaks = seq_len(simulation$settings$k)) +
    ggplot2::labs(x = "Grid index", y = "Centered feature effect",
                  color = NULL, linetype = NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    theme_pdp() +
    ggplot2::theme(legend.position = "bottom")

  masses <- simulation$observed_group_mass
  p_mass <- ggplot2::ggplot(masses, ggplot2::aes(group, mass)) +
    ggplot2::geom_hline(yintercept = simulation$settings$epsilon,
                        linetype = "dashed", color = pdp_palette[["blue"]]) +
    ggplot2::geom_col(fill = pdp_palette[["mid_grey"]], width = 0.72) +
    ggplot2::annotate(
      "text", x = 1, y = simulation$settings$epsilon,
      label = expression(epsilon), hjust = 0, vjust = -0.5, size = 3
    ) +
    ggplot2::scale_x_continuous(breaks = seq_len(simulation$settings$k)) +
    ggplot2::labs(x = "Context group", y = "Observed mass") +
    theme_pdp()

  summary <- summarize_s24(simulation)
  diagnostic_order <- c("Original P", "Average constraint",
                        "Pointwise constraints", "Intersection constraint")
  diagnostic_colors <- c(
    `Original P` = pdp_palette[["mid_grey"]],
    `Average constraint` = pdp_palette[["gold"]],
    `Pointwise constraints` = pdp_palette[["blue"]],
    `Intersection constraint` = pdp_palette[["purple"]]
  )
  coverage <- summary$diagnostics
  coverage$reference <- factor(coverage$reference, levels = diagnostic_order)
  coverage_long <- rbind(
    data.frame(reference = coverage$reference, value = coverage$min_coverage,
               diagnostic = "Minimum"),
    data.frame(reference = coverage$reference, value = coverage$mean_coverage,
               diagnostic = "Mean"),
    data.frame(reference = coverage$reference,
               value = coverage$simultaneous_valid_mass,
               diagnostic = "Simultaneous")
  )
  coverage_long <- coverage_long[is.finite(coverage_long$value), , drop = FALSE]
  zero_simultaneous <- coverage_long[
    coverage_long$diagnostic == "Simultaneous" &
      coverage_long$value == 0, , drop = FALSE
  ]
  p_coverage <- ggplot2::ggplot(
    coverage_long, ggplot2::aes(reference, value, fill = diagnostic)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78),
                      width = 0.7) +
    ggplot2::geom_point(
      data = zero_simultaneous,
      ggplot2::aes(reference, value), inherit.aes = FALSE,
      shape = 23, size = 2.6, stroke = 0.7,
      color = pdp_palette[["purple"]], fill = "white"
    ) +
    ggplot2::geom_text(
      data = zero_simultaneous,
      ggplot2::aes(reference, value, label = "0"), inherit.aes = FALSE,
      nudge_y = 0.035, size = 2.7, color = pdp_palette[["purple"]]
    ) +
    ggplot2::annotate(
      "text", x = 4, y = 0.055, label = "infeasible", size = 2.6,
      color = pdp_palette[["purple"]]
    ) +
    ggplot2::scale_fill_manual(values = c(
      Minimum = pdp_palette[["blue"]], Mean = pdp_palette[["gold"]],
      Simultaneous = pdp_palette[["purple"]]
    )) +
    ggplot2::scale_x_discrete(drop = FALSE, labels = c(
      `Original P` = "Original P",
      `Average constraint` = "Average-only",
      `Pointwise constraints` = "RCPD",
      `Intersection constraint` = "Intersection"
    )) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "Validity mass", fill = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    theme_pdp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
                   legend.position = "bottom")

  errors <- summary$metrics
  errors$method <- factor(errors$method, levels = method_order[-1])
  p_error <- ggplot2::ggplot(
    errors, ggplot2::aes(method, normalized_shape_error, color = method)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, normalized_shape_error - 2 * normalized_shape_error_mc_se),
                   ymax = normalized_shape_error + 2 * normalized_shape_error_mc_se),
      width = 0.15, linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_discrete(labels = c(
      PDP = "PDP", `Pointwise trim` = "PTPD",
      `Pointwise constraints` = "RCPD",
      `Average constraint` = "Average-only"
    )) +
    ggplot2::labs(x = NULL, y = "Normalized shape error") +
    theme_pdp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                   legend.position = "none")

  paths <- file.path(
    directory,
    sprintf("s24-constraint-comparison-%s.png", letters[1:4])
  )
  save_pdp_plot(p_curves, paths[[1]], width = 5.2, height = 4.0)
  save_pdp_plot(p_mass, paths[[2]], width = 5.2, height = 4.0)
  save_pdp_plot(p_coverage, paths[[3]], width = 5.2, height = 4.0)
  save_pdp_plot(p_error, paths[[4]], width = 5.2, height = 4.0)
  paths
}
