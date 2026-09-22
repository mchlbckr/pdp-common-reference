# Simulation S30: finite-sample Lipschitz sensitivity envelopes for the smooth
# interaction GAM in the held-out wine application.

s30_empirical_lipschitz_constant <- function(predictions, contexts,
                                              tolerance = 1e-10) {
  distances <- euclidean_cross_distance(contexts, contexts)
  gaps <- abs(outer(predictions, predictions, "-"))
  duplicate_conflict <- distances <= tolerance & gaps > tolerance
  if (any(duplicate_conflict)) {
    stop("duplicate contexts have conflicting predictions.", call. = FALSE)
  }
  ratios <- gaps[distances > tolerance] / distances[distances > tolerance]
  if (length(ratios) == 0L) 0 else max(ratios)
}

run_s30 <- function(case_study, model_name = "interaction_gam",
                    multipliers = c(1, 1.5, 2),
                    prediction_bounds = c(0, 10),
                    context_names = c(
                      "volatile.acidity", "sulphates", "density",
                      "total.sulfur.dioxide", "chlorides"
                    )) {
  validity <- attr(case_study, "evaluation_validity")
  prediction_matrices <- attr(case_study, "all_prediction_matrices")
  if (is.null(prediction_matrices)) {
    prediction_matrices <- attr(case_study, "prediction_matrices")
  }
  evaluation_contexts <- attr(case_study, "evaluation_contexts")
  if (is.null(validity) || is.null(prediction_matrices) ||
      is.null(evaluation_contexts)) {
    stop("the wine result lacks held-out matrices or contexts.", call. = FALSE)
  }
  predictions <- prediction_matrices[[model_name]]
  if (is.null(predictions)) stop("unknown wine model name.", call. = FALSE)
  grid <- sort(unique(case_study$z))
  if (!all(context_names %in% names(evaluation_contexts))) {
    stop("requested context variables are absent from the wine data.", call. = FALSE)
  }
  raw_contexts <- evaluation_contexts[, context_names, drop = FALSE]
  scales <- vapply(raw_contexts, stats::sd, numeric(1))
  keep <- is.finite(scales) & scales > 0
  standardized <- scale(
    raw_contexts[, keep, drop = FALSE], center = TRUE, scale = scales[keep]
  )

  pointwise_constants <- vapply(seq_along(grid), function(index) {
    selected <- validity[, index]
    s30_empirical_lipschitz_constant(
      predictions[selected, index], standardized[selected, , drop = FALSE]
    )
  }, numeric(1))
  admissible_minimum <- max(pointwise_constants)
  lipschitz_values <- admissible_minimum * multipliers
  n <- nrow(validity)

  range_rows <- lapply(seq_along(grid), function(index) {
    selected <- validity[, index]
    known <- sum(predictions[selected, index]) / n
    invalid_mass <- 1 - mean(selected)
    data.frame(
      z = grid[[index]], scenario = "Range only",
      lower = known + prediction_bounds[[1]] * invalid_mass,
      upper = known + prediction_bounds[[2]] * invalid_mass,
      isc = mean(selected), multiplier = Inf
    )
  })
  lipschitz_rows <- lapply(seq_along(grid), function(index) {
    selected <- validity[, index]
    valid_predictions <- predictions[selected, index]
    distances <- euclidean_cross_distance(
      standardized, standardized[selected, , drop = FALSE]
    )
    do.call(rbind, lapply(seq_along(lipschitz_values), function(l_index) {
      constant <- lipschitz_values[[l_index]]
      lower_candidates <- sweep(-constant * distances, 2,
                                valid_predictions, "+")
      upper_candidates <- sweep(constant * distances, 2,
                                valid_predictions, "+")
      lower_envelope <- pmax(
        prediction_bounds[[1]], apply(lower_candidates, 1, max)
      )
      upper_envelope <- pmin(
        prediction_bounds[[2]], apply(upper_candidates, 1, min)
      )
      data.frame(
        z = grid[[index]],
        scenario = sprintf("L = %.1f L0", multipliers[[l_index]]),
        lower = mean(lower_envelope), upper = mean(upper_envelope),
        isc = mean(selected), multiplier = multipliers[[l_index]]
      )
    }))
  })
  bounds <- rbind(do.call(rbind, lipschitz_rows), do.call(rbind, range_rows))
  bounds$width <- bounds$upper - bounds$lower
  bounds$midpoint <- (bounds$lower + bounds$upper) / 2
  bounds$pdp <- colMeans(predictions)[match(bounds$z, grid)]
  scenario_order <- c(
    sprintf("L = %.1f L0", multipliers), "Range only"
  )
  bounds$scenario <- factor(bounds$scenario, levels = scenario_order)
  list(
    bounds = bounds, pointwise_constants = pointwise_constants,
    admissible_minimum = admissible_minimum, multipliers = multipliers,
    prediction_bounds = prediction_bounds, context_scales = scales[keep],
    metric = "Euclidean distance after standardizing each context variable by its held-out standard deviation",
    model_name = model_name,
    interpretation = paste(
      "finite-sample sensitivity envelope for a smooth fitted GAM;",
      "the empirical anchor does not establish a population Lipschitz constant"
    )
  )
}

summarize_s30 <- function(simulation) {
  groups <- split(simulation$bounds, simulation$bounds$scenario)
  widths <- do.call(rbind, lapply(groups, function(part) {
    data.frame(
      scenario = as.character(part$scenario[[1]]),
      mean_width = mean(part$width), max_width = max(part$width),
      min_isc = min(part$isc), stringsAsFactors = FALSE
    )
  }))
  lf0 <- simulation$bounds[
    as.character(simulation$bounds$scenario) == "L = 1.0 L0", , drop = FALSE
  ]
  lower_gap <- lf0$pdp - lf0$lower
  upper_gap <- lf0$upper - lf0$pdp
  list(
    admissible_minimum = simulation$admissible_minimum,
    widths = widths, metric = simulation$metric,
    prediction_bounds = simulation$prediction_bounds,
    model_name = simulation$model_name,
    interpretation = simulation$interpretation,
    pdp_inside_lf0 = all(lower_gap >= -1e-10 & upper_gap >= -1e-10),
    minimum_pdp_band_margin = min(lower_gap, upper_gap),
    maximum_pdp_band_violation = max(pmax(-lower_gap, -upper_gap, 0))
  )
}

plot_s30 <- function(simulation,
                     path = file.path("results", "s30-lipschitz-wine.png")) {
  bounds <- simulation$bounds
  scenario_order <- levels(bounds$scenario)
  colors <- c(
    `L = 1.0 L0` = pdp_palette[["blue"]],
    `L = 1.5 L0` = pdp_palette[["orange"]],
    `L = 2.0 L0` = pdp_palette[["purple"]],
    `Range only` = pdp_palette[["mid_grey"]]
  )
  facet_labels <- c(
    `L = 1.0 L0` = "L[f] == 1.0*L[f*','*0]",
    `L = 1.5 L0` = "L[f] == 1.5*L[f*','*0]",
    `L = 2.0 L0` = "L[f] == 2.0*L[f*','*0]",
    `Range only` = "Range~only"
  )
  legend_labels <- c(
    expression(L[f] == 1.0 * L[f,0]),
    expression(L[f] == 1.5 * L[f,0]),
    expression(L[f] == 2.0 * L[f,0]),
    expression("Range only")
  )
  p_bands <- ggplot2::ggplot(
    bounds, ggplot2::aes(z, ymin = lower, ymax = upper, fill = scenario)
  ) +
    ggplot2::geom_ribbon(alpha = 0.24, color = NA) +
    ggplot2::geom_line(
      ggplot2::aes(y = lower, color = scenario), linewidth = 0.45
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = upper, color = scenario), linewidth = 0.45
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = pdp), color = pdp_palette[["ink"]],
      linetype = "dashed", linewidth = 0.65
    ) +
    ggplot2::facet_wrap(
      ~scenario, ncol = 2,
      labeller = ggplot2::labeller(
        scenario = ggplot2::as_labeller(facet_labels, ggplot2::label_parsed)
      )
    ) +
    ggplot2::scale_fill_manual(values = colors, labels = legend_labels,
                               drop = FALSE) +
    ggplot2::scale_color_manual(values = colors, labels = legend_labels,
                                drop = FALSE) +
    ggplot2::labs(x = "Alcohol (% by volume)",
                  y = "Predicted wine quality") +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "none")

  p_width <- ggplot2::ggplot(
    bounds, ggplot2::aes(z, width, color = scenario, linetype = scenario)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors, labels = legend_labels,
                                drop = FALSE) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(c("solid", "dashed", "dotdash", "dotted"),
                               scenario_order), labels = legend_labels,
      drop = FALSE
    ) +
    ggplot2::labs(x = "Alcohol (% by volume)", y = "Robustness-envelope width",
                  color = NULL, linetype = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  save_pdp_panel_set(
    list(p_bands, p_width), path,
    widths = 10, heights = c(5.3, 3.6)
  )
}
