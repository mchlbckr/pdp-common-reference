# Simulation S29: an externally anchored regime in which relaxation is the
# most accurate support-qualified fixed-reference compromise.

s29_center <- function(x) x - mean(x)

s29_signal <- function(z, half_width = 0.5) {
  sin(pi * z / (2 * half_width))
}

simulate_s29_data <- function(n, context_mean = 0.2, context_sd = 0.7,
                              residual_sd = 0.35, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- stats::rnorm(n, context_mean, context_sd)
  data.frame(x = u + stats::rnorm(n, 0, residual_sd), u = u)
}

s29_validity <- function(residual_sd = 0.35, alpha = 0.05) {
  radius <- stats::qnorm(1 - alpha / 2) * residual_sd
  function(newdata) abs(newdata$x - newdata$u) <= radius
}

s29_prediction <- function(half_width = 0.5, beta = 1.4,
                           interaction = 2.5) {
  function(newdata) {
    s29_signal(newdata$x, half_width) + beta * newdata$u +
      interaction * newdata$x * newdata$u
  }
}

s29_one <- function(repetition, n_evaluation, context_mean, context_sd,
                    residual_sd, alpha, epsilon, grid, beta, interaction,
                    ale_bins, seed) {
  data <- simulate_s29_data(
    n_evaluation, context_mean, context_sd, residual_sd, seed + repetition
  )
  half_width <- max(abs(grid))
  is_valid <- s29_validity(residual_sd, alpha)
  predict_function <- s29_prediction(half_width, beta, interaction)
  validity <- validity_matrix(data, "x", grid, is_valid)
  common <- rowSums(validity) == ncol(validity)
  if (!any(common)) {
    stop("S29 requires at least one hard-common evaluation context.", call. = FALSE)
  }
  relaxed <- fit_relaxed_reference(validity, epsilon, tolerance = 2e-5)

  pd <- estimate_pdp(data, "x", grid, predict_function)$pdp
  ptpd <- estimate_pointwise_trimmed_pdp(
    data, "x", grid, is_valid, predict_function
  )$pointwise_trimmed_pdp
  hard <- estimate_weighted_pdp(
    data, "x", grid, predict_function, common / sum(common)
  )$weighted_pdp
  rcpd <- estimate_weighted_pdp(
    data, "x", grid, predict_function, relaxed$weights
  )$weighted_pdp
  breaks <- unique(stats::quantile(
    data$x, seq(0, 1, length.out = ale_bins + 1L), names = FALSE
  ))
  ale_raw <- estimate_ale(data, "x", breaks, predict_function)
  ale <- stats::approx(ale_raw$z, ale_raw$ale, xout = grid, rule = 2)$y

  target <- s29_signal(grid, half_width) + interaction * grid * context_mean
  target_centered <- s29_center(target)
  target_amplitude <- diff(range(target_centered))
  curves <- list(
    `Ordinary PD` = pd, PTPD = ptpd, `Hard CSPD` = hard,
    RCPD = rcpd, ALE = ale
  )
  metrics <- do.call(rbind, lapply(names(curves), function(method) {
    estimate <- s29_center(curves[[method]])
    data.frame(
      repetition = repetition, method = method,
      normalized_shape_error = max(abs(estimate - target_centered)) /
        target_amplitude,
      normalized_rmse = sqrt(mean((estimate - target_centered)^2)) /
        target_amplitude,
      stringsAsFactors = FALSE
    )
  }))
  curve_data <- do.call(rbind, lapply(names(curves), function(method) {
    data.frame(
      repetition = repetition, z = grid,
      value = s29_center(curves[[method]]), method = method,
      stringsAsFactors = FALSE
    )
  }))
  pointwise_context <- vapply(grid, function(z) {
    selected <- is_valid(make_queries(data, "x", z))
    mean(data$u[selected])
  }, numeric(1))
  context_data <- rbind(
    data.frame(repetition = repetition, z = grid, value = pointwise_context,
               reference = "PTPD"),
    data.frame(repetition = repetition, z = grid, value = mean(data$u),
               reference = "Original P"),
    data.frame(repetition = repetition, z = grid, value = mean(data$u[common]),
               reference = "Hard CSPD"),
    data.frame(repetition = repetition, z = grid,
               value = sum(relaxed$weights * data$u), reference = "RCPD")
  )
  diagnostics <- data.frame(
    repetition = repetition, hard_common_mass = mean(common),
    original_min_coverage = min(colMeans(validity)),
    rcpd_min_coverage = relaxed$min_coverage,
    rcpd_ess_fraction = relaxed$effective_sample_size / n_evaluation,
    rcpd_simultaneous_mass = relaxed$simultaneous_valid_mass,
    stringsAsFactors = FALSE
  )
  list(metrics = metrics, curves = curve_data, contexts = context_data,
       diagnostics = diagnostics,
       target = data.frame(z = grid, value = target_centered))
}

run_s29 <- function(repetitions = 300L, n_evaluation = 400L,
                    context_mean = 0.2, context_sd = 0.7,
                    residual_sd = 0.35, alpha = 0.05, epsilon = 0.30,
                    grid = seq(-0.5, 0.5, length.out = 21L), beta = 1.4,
                    interaction = 2.5, ale_bins = 12L, seed = 20261001L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  grid <- assert_grid(grid)
  evaluate <- function(index) s29_one(
    index, n_evaluation, context_mean, context_sd, residual_sd, alpha,
    epsilon, grid, beta, interaction, ale_bins, seed
  )
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(repetitions), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(repetitions), evaluate)
  }
  list(
    metrics = do.call(rbind, lapply(pieces, `[[`, "metrics")),
    curves = do.call(rbind, lapply(pieces, `[[`, "curves")),
    contexts = do.call(rbind, lapply(pieces, `[[`, "contexts")),
    diagnostics = do.call(rbind, lapply(pieces, `[[`, "diagnostics")),
    target = pieces[[1]]$target,
    settings = list(
      repetitions = repetitions, n_evaluation = n_evaluation,
      context_mean = context_mean, context_sd = context_sd,
      residual_sd = residual_sd, alpha = alpha, epsilon = epsilon,
      grid = grid, beta = beta, interaction = interaction,
      ale_bins = ale_bins
    )
  )
}

summarize_s29 <- function(simulation) {
  metric_groups <- split(simulation$metrics, simulation$metrics$method)
  metrics <- do.call(rbind, lapply(metric_groups, function(part) {
    data.frame(
      method = part$method[[1]],
      normalized_shape_error = mean(part$normalized_shape_error),
      normalized_shape_error_mc_se = stats::sd(part$normalized_shape_error) /
        sqrt(nrow(part)),
      normalized_rmse = mean(part$normalized_rmse),
      stringsAsFactors = FALSE
    )
  }))
  curve_groups <- split(
    simulation$curves,
    interaction(simulation$curves$method, simulation$curves$z, drop = TRUE)
  )
  curves <- do.call(rbind, lapply(curve_groups, function(part) {
    data.frame(
      method = part$method[[1]], z = part$z[[1]], value = mean(part$value),
      mc_se = stats::sd(part$value) / sqrt(nrow(part)),
      stringsAsFactors = FALSE
    )
  }))
  context_groups <- split(
    simulation$contexts,
    interaction(simulation$contexts$reference, simulation$contexts$z,
                drop = TRUE)
  )
  contexts <- do.call(rbind, lapply(context_groups, function(part) {
    data.frame(
      reference = part$reference[[1]], z = part$z[[1]],
      value = mean(part$value),
      mc_se = stats::sd(part$value) / sqrt(nrow(part)),
      stringsAsFactors = FALSE
    )
  }))
  diagnostic_summary <- data.frame(
    hard_common_mass = mean(simulation$diagnostics$hard_common_mass),
    original_min_coverage = mean(simulation$diagnostics$original_min_coverage),
    rcpd_min_coverage = mean(simulation$diagnostics$rcpd_min_coverage),
    rcpd_ess_fraction = mean(simulation$diagnostics$rcpd_ess_fraction),
    rcpd_simultaneous_mass = mean(
      simulation$diagnostics$rcpd_simultaneous_mass
    )
  )
  list(metrics = metrics, curves = curves, contexts = contexts,
       diagnostics = diagnostic_summary, settings = simulation$settings,
       target = simulation$target)
}

plot_s29 <- function(summary,
                     path = file.path("results", "s29-rcpd-comparison.png")) {
  method_order <- c("Target", "Ordinary PD", "PTPD", "Hard CSPD", "RCPD", "ALE")
  colors <- c(
    Target = pdp_palette[["ink"]], `Ordinary PD` = pdp_palette[["mid_grey"]],
    PTPD = pdp_palette[["purple"]], `Hard CSPD` = pdp_palette[["green"]],
    RCPD = pdp_palette[["blue"]], ALE = pdp_palette[["gold"]]
  )
  line_types <- c(
    Target = "solid", `Ordinary PD` = "dashed", PTPD = "dotdash",
    `Hard CSPD` = "longdash", RCPD = "solid", ALE = "dotted"
  )
  curves <- rbind(
    transform(summary$target, method = "Target", mc_se = 0),
    summary$curves
  )
  curves$method <- factor(curves$method, levels = method_order)
  p_curves <- ggplot2::ggplot(
    curves, ggplot2::aes(z, value, color = method, linetype = method)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = line_types) +
    ggplot2::labs(x = expression(z), y = "Centered feature effect",
                  color = NULL, linetype = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  context_order <- c("Original P", "PTPD", "Hard CSPD", "RCPD")
  context_colors <- c(
    `Original P` = pdp_palette[["mid_grey"]],
    PTPD = pdp_palette[["purple"]],
    `Hard CSPD` = pdp_palette[["green"]],
    RCPD = pdp_palette[["blue"]]
  )
  context_line_types <- c(
    `Original P` = "dashed", PTPD = "dotdash",
    `Hard CSPD` = "longdash", RCPD = "solid"
  )
  contexts <- summary$contexts
  contexts$reference <- factor(contexts$reference, levels = context_order)
  p_context <- ggplot2::ggplot(
    contexts, ggplot2::aes(z, value, color = reference, linetype = reference)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(
      values = context_colors, breaks = context_order, labels = context_order
    ) +
    ggplot2::scale_linetype_manual(
      values = context_line_types, breaks = context_order, labels = context_order
    ) +
    ggplot2::labs(x = expression(z), y = "Mean context U",
                  color = NULL, linetype = NULL) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    theme_pdp(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")

  errors <- summary$metrics
  errors$method <- factor(errors$method, levels = method_order[-1])
  p_error <- ggplot2::ggplot(
    errors, ggplot2::aes(method, normalized_shape_error, color = method)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, normalized_shape_error - 2 * normalized_shape_error_mc_se),
        ymax = normalized_shape_error + 2 * normalized_shape_error_mc_se
      ), width = 0.14, linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", normalized_shape_error)),
      nudge_y = 0.015, size = 2.7, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = NULL, y = "Normalized shape error") +
    theme_pdp(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
                   legend.position = "none")

  save_pdp_panel_set(
    list(p_curves, p_context, p_error), path,
    widths = 5.2, heights = 4.0
  )
}
