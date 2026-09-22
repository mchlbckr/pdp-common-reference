# Simulation S18: hard versus relaxed common references.

s18_prediction <- function(newdata) {
  sin(newdata$x1) + 0.5 * newdata$x2^2 + 0.9 * newdata$x1 * newdata$x2 +
    1.2 * newdata$x1^2 * newdata$x2^2
}

s18_one <- function(repetition, n_evaluation, rho, grid_half_width, grid_points,
                    alpha, epsilons, seed) {
  set.seed(seed + repetition)
  evaluation <- simulate_s1_data(n_evaluation, rho)
  grid <- seq(-grid_half_width, grid_half_width, length.out = grid_points)
  valid <- s14_oracle_validity(rho, alpha)
  validity <- validity_matrix(evaluation, "x1", grid, valid)
  pdp <- estimate_pdp(evaluation, "x1", grid, s18_prediction)$pdp
  hard <- estimate_relaxed_cspd(evaluation, "x1", grid, valid, s18_prediction, epsilon = 0)

  rows <- do.call(rbind, lapply(epsilons, function(epsilon) {
    relaxed <- estimate_relaxed_cspd(
      evaluation, "x1", grid, valid, s18_prediction,
      epsilon = epsilon, tolerance = 2e-5
    )
    data.frame(
      repetition = repetition, epsilon = epsilon,
      hard_gamma = unique(hard$effective_sample_size) / n_evaluation,
      mean_isc = mean(validity),
      observational_acceptance = mean(valid(evaluation)),
      min_coverage = unique(relaxed$min_coverage),
      mean_coverage = unique(relaxed$mean_coverage),
      simultaneous_valid_mass = unique(relaxed$simultaneous_valid_mass),
      kl = unique(relaxed$kl),
      effective_sample_size = unique(relaxed$effective_sample_size),
      effective_sample_fraction = unique(relaxed$effective_sample_size) / n_evaluation,
      observed_patterns = unique(relaxed$observed_patterns),
      median_pattern_size = unique(relaxed$median_pattern_size),
      min_pattern_size = unique(relaxed$min_pattern_size),
      singleton_pattern_fraction = unique(relaxed$singleton_pattern_fraction),
      active_constraints = unique(relaxed$active_constraints),
      active_affine_dimension = unique(relaxed$active_affine_dimension),
      active_licq = unique(relaxed$active_licq),
      binding_constraints = unique(relaxed$binding_constraints),
      strict_complementarity = unique(relaxed$strict_complementarity),
      max_weight = unique(relaxed$max_weight),
      epsilon_min = unique(relaxed$epsilon_min),
      max_abs_weight_prediction_correlation = max(abs(relaxed$weight_prediction_correlation)),
      max_covariance_identity_error = max(relaxed$covariance_identity_error),
      covariance_path_range = diff(range(relaxed$reference_covariance)),
      centered_distance_to_pdp = centered_curve_distance(relaxed$relaxed_cspd, pdp),
      centered_distance_to_hard = centered_curve_distance(relaxed$relaxed_cspd, hard$relaxed_cspd),
      stringsAsFactors = FALSE
    )
  }))

  curves <- if (repetition == 1L) {
    relaxed_curves <- do.call(rbind, lapply(epsilons, function(epsilon) {
      estimate <- estimate_relaxed_cspd(
        evaluation, "x1", grid, valid, s18_prediction,
        epsilon = epsilon, tolerance = 2e-5
      )
      data.frame(z = grid, value = estimate$relaxed_cspd,
                 curve = if (epsilon == 0) "Hard CSPD" else sprintf("Relaxed, epsilon = %.2f", epsilon))
    }))
    rbind(data.frame(z = grid, value = pdp, curve = "PDP"), relaxed_curves)
  } else NULL

  list(metrics = rows, curves = curves,
       sample_focal = if (repetition == 1L) evaluation$x1 else NULL)
}

run_s18 <- function(repetitions = 100L, n_evaluation = 1000L, rho = 0.8,
                    grid_half_width = 1.1, grid_points = 21L, alpha = 0.05,
                    epsilons = c(0, 0.05, 0.1, 0.2, 0.3), seed = 20260923L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  indices <- seq_len(repetitions)
  evaluate <- function(index) s18_one(index, n_evaluation, rho, grid_half_width,
                                      grid_points, alpha, epsilons, seed)
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(indices, evaluate, mc.cores = workers)
  } else {
    lapply(indices, evaluate)
  }
  list(
    metrics = do.call(rbind, lapply(pieces, `[[`, "metrics")),
    curves = pieces[[1]]$curves,
    sample_focal = pieces[[1]]$sample_focal,
    settings = list(repetitions = repetitions, n_evaluation = n_evaluation,
                    rho = rho, grid_half_width = grid_half_width,
                    grid_points = grid_points, alpha = alpha, seed = seed)
  )
}

summarize_s18 <- function(simulation) {
  measures <- c("hard_gamma", "mean_isc", "observational_acceptance", "min_coverage",
                "mean_coverage", "simultaneous_valid_mass",
                "kl", "effective_sample_size", "effective_sample_fraction",
                "centered_distance_to_pdp", "centered_distance_to_hard",
                "observed_patterns", "median_pattern_size", "min_pattern_size",
                "singleton_pattern_fraction", "active_constraints",
                "active_affine_dimension", "active_licq",
                "binding_constraints", "strict_complementarity", "max_weight",
                "epsilon_min", "max_abs_weight_prediction_correlation",
                "max_covariance_identity_error", "covariance_path_range")
  pieces <- split(simulation$metrics, simulation$metrics$epsilon)
  do.call(rbind, lapply(pieces, function(part) {
    output <- data.frame(epsilon = part$epsilon[[1]])
    for (measure in measures) {
      output[[measure]] <- mean(part[[measure]])
      output[[paste0(measure, "_mc_se")]] <- stats::sd(part[[measure]]) / sqrt(nrow(part))
    }
    output
  }))
}

plot_s18 <- function(simulation, path = file.path("results", "s18-relaxed-reference.png")) {
  summary <- summarize_s18(simulation)
  p_kl <- ggplot2::ggplot(summary, ggplot2::aes(epsilon, kl)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = kl - 2 * kl_mc_se, ymax = kl + 2 * kl_mc_se),
      width = 0.012, linewidth = 0.6, color = pdp_palette[["blue"]]
    ) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.7) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.2) +
    ggplot2::labs(x = expression(epsilon), y = "KL distortion") +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  p_ess <- ggplot2::ggplot(summary, ggplot2::aes(epsilon, effective_sample_fraction)) +
    ggplot2::geom_hline(yintercept = summary$hard_gamma[[1]], linetype = "dashed",
                        color = pdp_palette[["mid_grey"]], linewidth = 0.55) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.7) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.2) +
    ggplot2::annotate("text", x = max(summary$epsilon), y = summary$hard_gamma[[1]],
                      label = "Hard common mass", hjust = 1, vjust = -0.6,
                      size = 2.8, color = pdp_palette[["mid_grey"]]) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = expression(epsilon), y = "Effective sample fraction") +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  selected <- simulation$curves$curve %in% c("PDP", "Hard CSPD", "Relaxed, epsilon = 0.10", "Relaxed, epsilon = 0.20")
  curves <- simulation$curves[selected, ]
  curves$value <- ave(curves$value, curves$curve, FUN = function(x) x - mean(x))
  curve_order <- c("PDP", "Hard CSPD", "Relaxed, epsilon = 0.10",
                   "Relaxed, epsilon = 0.20")
  curve_labels <- c("PDP", "Hard CSPD", expression(RCPD~(epsilon==0.10)),
                    expression(RCPD~(epsilon==0.20)))
  curves$curve <- factor(curves$curve, levels = curve_order)
  rug <- data.frame(z = simulation$sample_focal)
  rug <- rug[rug$z >= min(curves$z) & rug$z <= max(curves$z), , drop = FALSE]
  curve_colors <- c("PDP" = pdp_palette[["ink"]],
                    "Hard CSPD" = pdp_palette[["orange"]],
                    "Relaxed, epsilon = 0.10" = pdp_palette[["blue"]],
                    "Relaxed, epsilon = 0.20" = pdp_palette[["green"]])
  p_curves <- ggplot2::ggplot(curves, ggplot2::aes(z, value, color = curve,
                                                   linetype = curve)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_rug(data = rug,
                      ggplot2::aes(x = z), inherit.aes = FALSE, sides = "b",
                      alpha = 0.14, color = pdp_palette[["mid_grey"]]) +
    ggplot2::scale_color_manual(values = curve_colors, labels = curve_labels) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted"),
                                   labels = curve_labels) +
    ggplot2::labs(x = expression(z), y = "Centered feature-effect curve",
                  color = NULL, linetype = NULL) +
    theme_pdp()
  p_distance <- ggplot2::ggplot(summary, ggplot2::aes(epsilon, centered_distance_to_pdp)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = centered_distance_to_pdp - 2 * centered_distance_to_pdp_mc_se,
                   ymax = centered_distance_to_pdp + 2 * centered_distance_to_pdp_mc_se),
      width = 0.012, linewidth = 0.6, color = pdp_palette[["blue"]]
    ) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.7) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.2) +
    ggplot2::labs(x = expression(epsilon), y = "Centered distance to PDP") +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  save_pdp_panel_set(
    list(p_kl, p_ess, p_curves, p_distance), path,
    widths = 5.2, heights = 4.0
  )
}
