# Simulation S32: a fixed off-support perturbation on real Bank
# Marketing contexts. The perturbation is semisynthetic and changes only
# queries rejected by the fitted validity rule.

s32_scale_focal <- function(z, grid) {
  midpoint <- mean(range(grid))
  half_width <- diff(range(grid)) / 2
  if (!is.finite(half_width) || half_width <= 0) {
    stop("grid must span a nonzero interval.", call. = FALSE)
  }
  pmax(-1, pmin(1, (z - midpoint) / half_width))
}

s32_perturb_probabilities <- function(predictions, invalid, focal_scale,
                                      logit_shift) {
  if (!is.numeric(predictions) || any(!is.finite(predictions)) ||
      any(predictions < 0 | predictions > 1) ||
      !is.logical(invalid) || length(invalid) != length(predictions) ||
      length(focal_scale) != length(predictions)) {
    stop("predictions, invalid flags, and focal scales are incompatible.",
         call. = FALSE)
  }
  output <- predictions
  if (any(invalid)) {
    clipped <- pmin(pmax(predictions[invalid], 1e-6), 1 - 1e-6)
    output[invalid] <- stats::plogis(
      stats::qlogis(clipped) + logit_shift * focal_scale[invalid]
    )
  }
  output
}

run_s32 <- function(bank_case, logit_shift = 3, epsilon = 0.10) {
  predictions <- attr(bank_case, "prediction_matrix")
  validity <- attr(bank_case, "evaluation_validity")
  observed_predictions <- attr(bank_case, "observed_predictions")
  observed_validity <- attr(bank_case, "observed_validity")
  observed_response <- attr(bank_case, "evaluation_response")
  observed_focal <- attr(bank_case, "sample_focal")
  grid <- sort(unique(bank_case$z))
  if (!is.matrix(predictions) || !is.matrix(validity) ||
      !identical(dim(predictions), dim(validity)) ||
      ncol(predictions) != length(grid) || !is.logical(validity)) {
    stop("bank case lacks compatible prediction and validity matrices.",
         call. = FALSE)
  }
  if (!is.numeric(logit_shift) || length(logit_shift) != 1L ||
      !is.finite(logit_shift) || logit_shift <= 0) {
    stop("logit_shift must be one positive finite number.", call. = FALSE)
  }
  grid_scale <- s32_scale_focal(grid, grid)
  shift_scale <- matrix(
    rep(grid_scale, each = nrow(predictions)), nrow = nrow(predictions)
  )
  stressed <- s32_perturb_probabilities(
    predictions, !validity, shift_scale, logit_shift
  )

  relaxed <- fit_relaxed_reference(validity, epsilon = epsilon,
                                    tolerance = 2e-5)
  common <- apply(validity, 1, all)
  if (!any(common)) stop("stress test requires hard-common contexts.", call. = FALSE)
  pointwise_mean <- function(matrix, mask) {
    vapply(seq_len(ncol(matrix)), function(index) {
      mean(matrix[mask[, index], index])
    }, numeric(1))
  }
  curves_for <- function(matrix) {
    list(
      `Ordinary PD` = colMeans(matrix),
      RCPD = colSums(matrix * relaxed$weights),
      `Hard CSPD` = colMeans(matrix[common, , drop = FALSE]),
      PTPD = pointwise_mean(matrix, validity)
    )
  }
  original_curves <- curves_for(predictions)
  stressed_curves <- curves_for(stressed)
  methods <- names(original_curves)
  center <- function(x) x - mean(x)
  curves <- do.call(rbind, lapply(methods, function(method) {
    rbind(
      data.frame(z = grid, value = center(original_curves[[method]]),
                 method = method, fit = "Original fit"),
      data.frame(z = grid, value = center(stressed_curves[[method]]),
                 method = method, fit = "Off-support perturbation")
    )
  }))
  sensitivity <- data.frame(
    method = methods,
    centered_sup_distance = vapply(methods, function(method) {
      max(abs(center(original_curves[[method]]) -
                center(stressed_curves[[method]])))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  observed_stressed <- NULL
  unchanged_fraction <- NA_real_
  original_brier <- NA_real_
  stressed_brier <- NA_real_
  if (is.numeric(observed_predictions) && is.logical(observed_validity) &&
      (is.numeric(observed_response) || is.logical(observed_response)) &&
      is.numeric(observed_focal) &&
      length(observed_predictions) == length(observed_validity) &&
      length(observed_predictions) == length(observed_response) &&
      length(observed_predictions) == length(observed_focal)) {
    observed_stressed <- s32_perturb_probabilities(
      observed_predictions, !observed_validity,
      s32_scale_focal(observed_focal, grid), logit_shift
    )
    unchanged_fraction <- mean(observed_stressed == observed_predictions)
    observed_response <- as.numeric(observed_response)
    original_brier <- mean((observed_response - observed_predictions)^2)
    stressed_brier <- mean((observed_response - observed_stressed)^2)
  }

  list(
    curves = curves, sensitivity = sensitivity,
    logit_shift = logit_shift, epsilon = epsilon,
    grid_invalid_fraction = mean(!validity),
    observed_unchanged_fraction = unchanged_fraction,
    original_brier = original_brier, stressed_brier = stressed_brier,
    brier_change = stressed_brier - original_brier,
    relaxed_min_coverage = relaxed$min_coverage
  )
}

summarize_s32 <- function(simulation) {
  list(
    sensitivity = simulation$sensitivity,
    logit_shift = simulation$logit_shift,
    epsilon = simulation$epsilon,
    grid_invalid_fraction = simulation$grid_invalid_fraction,
    observed_unchanged_fraction = simulation$observed_unchanged_fraction,
    original_brier = simulation$original_brier,
    stressed_brier = simulation$stressed_brier,
    brier_change = simulation$brier_change,
    relaxed_min_coverage = simulation$relaxed_min_coverage
  )
}

plot_s32 <- function(simulation,
                     path = file.path("results", "s32-bank-offsupport-stress.png")) {
  curves <- simulation$curves
  curves$method <- factor(
    curves$method, levels = c("Ordinary PD", "RCPD", "Hard CSPD", "PTPD")
  )
  fit_colors <- c(
    `Original fit` = pdp_palette[["ink"]],
    `Off-support perturbation` = pdp_palette[["purple"]]
  )
  p_curves <- ggplot2::ggplot(
    curves, ggplot2::aes(z, value, color = fit, linetype = fit)
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::facet_wrap(~method, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = fit_colors) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
    ggplot2::labs(x = "Age", y = "Centered predicted probability",
                  color = NULL, linetype = NULL) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  distances <- simulation$sensitivity
  distances$method <- factor(
    distances$method, levels = c("Ordinary PD", "RCPD", "Hard CSPD", "PTPD")
  )
  p_distance <- ggplot2::ggplot(
    distances, ggplot2::aes(method, centered_sup_distance, fill = method)
  ) +
    ggplot2::geom_col(width = 0.68, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", centered_sup_distance)),
      vjust = -0.4, size = 3
    ) +
    ggplot2::scale_fill_manual(values = c(
      `Ordinary PD` = pdp_palette[["purple"]], RCPD = pdp_palette[["blue"]],
      `Hard CSPD` = pdp_palette[["orange"]], PTPD = pdp_palette[["green"]]
    )) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(x = NULL, y = "Centered sensitivity") +
    theme_pdp(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  save_pdp_panel_set(
    list(p_curves, p_distance), path,
    widths = c(7, 4.3), heights = 4.5
  )
}
