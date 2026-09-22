# Simulation S23: pointwise trimming creates composition bias even when the
# supported prediction surface is additive.

s23_truth <- function(z) sin(pi * z / 2)

simulate_s23_data <- function(n, sigma = 0.45, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- stats::rnorm(n)
  data.frame(x = u + sigma * stats::rnorm(n), u = u)
}

s23_validity <- function(sigma = 0.45, alpha = 0.05) {
  radius <- stats::qnorm(1 - alpha / 2) * sigma
  function(newdata) abs(newdata$x - newdata$u) <= radius
}

s23_prediction <- function(direction = 0, sigma = 0.45, alpha = 0.05,
                           beta = 1.5, strength = 4) {
  radius <- stats::qnorm(1 - alpha / 2) * sigma
  function(newdata) {
    excess <- pmax(abs(newdata$x - newdata$u) - radius, 0)
    s23_truth(newdata$x) + beta * newdata$u +
      direction * strength * excess^2
  }
}

s23_center <- function(x) x - mean(x)

s23_one <- function(repetition, n_evaluation, sigma, alpha, epsilon, grid,
                    beta, strength, seed) {
  data <- simulate_s23_data(n_evaluation, sigma, seed + repetition)
  valid <- s23_validity(sigma, alpha)
  h <- validity_matrix(data, "x", grid, valid)
  relaxed <- fit_relaxed_reference(h, epsilon = epsilon, tolerance = 2e-7)
  hard <- fit_relaxed_reference(h, epsilon = 0)
  truth <- s23_truth(grid)
  truth_centered <- s23_center(truth)
  truth_amplitude <- diff(range(truth_centered))

  extension_names <- c("Downward extension", "Upward extension")
  directions <- c(-1, 1)
  results <- lapply(seq_along(directions), function(index) {
    prediction <- s23_prediction(
      directions[[index]], sigma, alpha, beta, strength
    )
    pd <- estimate_pdp(data, "x", grid, prediction)$pdp
    trimmed <- estimate_pointwise_trimmed_pdp(
      data, "x", grid, valid, prediction
    )$pointwise_trimmed_pdp
    hard_curve <- estimate_weighted_pdp(
      data, "x", grid, prediction, hard$weights
    )$weighted_pdp
    relaxed_curve <- estimate_weighted_pdp(
      data, "x", grid, prediction, relaxed$weights
    )$weighted_pdp
    breaks <- unique(stats::quantile(
      data$x, probs = seq(0, 1, length.out = 22), names = FALSE
    ))
    ale <- estimate_ale(data, "x", breaks, prediction)
    ale_curve <- stats::approx(ale$z, ale$ale, xout = grid, rule = 2)$y
    curves <- list(
      PDP = pd,
      `Pointwise trim` = trimmed,
      `Hard CSPD` = hard_curve,
      RCPD = relaxed_curve,
      ALE = ale_curve
    )
    metrics <- do.call(rbind, lapply(names(curves), function(method) {
      centered <- s23_center(curves[[method]])
      data.frame(
        repetition = repetition,
        extension = extension_names[[index]],
        method = method,
        shape_error = max(abs(centered - truth_centered)),
        rmse = sqrt(mean((centered - truth_centered)^2)),
        normalized_shape_error = max(abs(centered - truth_centered)) /
          truth_amplitude,
        gamma = mean(rowSums(h) == ncol(h)),
        mean_isc = mean(h),
        rcpd_min_coverage = relaxed$min_coverage,
        rcpd_simultaneous_mass = relaxed$simultaneous_valid_mass,
        rcpd_ess_fraction = relaxed$effective_sample_size / n_evaluation,
        stringsAsFactors = FALSE
      )
    }))
    curve_data <- do.call(rbind, lapply(names(curves), function(method) {
      data.frame(
        z = grid, value = s23_center(curves[[method]]), method = method,
        extension = extension_names[[index]], stringsAsFactors = FALSE
      )
    }))
    list(metrics = metrics, curves = curve_data)
  })

  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  paired <- split(metrics, metrics$method)
  sensitivity <- do.call(rbind, lapply(paired, function(part) {
    curves <- lapply(results, function(result) {
      result$curves$value[result$curves$method == part$method[[1]]]
    })
    data.frame(
      repetition = repetition,
      method = part$method[[1]],
      extension_sensitivity = max(abs(curves[[1]] - curves[[2]])),
      stringsAsFactors = FALSE
    )
  }))
  metrics <- merge(metrics, sensitivity, by = c("repetition", "method"),
                   all.x = TRUE, sort = FALSE)

  pointwise_context_mean <- vapply(grid, function(z) {
    selected <- valid(make_queries(data, "x", z))
    mean(data$u[selected])
  }, numeric(1))
  composition <- data.frame(
    z = grid,
    pointwise_mean_u = pointwise_context_mean,
    original_mean_u = mean(data$u),
    hard_mean_u = sum(hard$weights * data$u),
    relaxed_mean_u = sum(relaxed$weights * data$u)
  )
  list(
    metrics = metrics,
    curves = if (repetition == 1L) do.call(rbind, lapply(results, `[[`, "curves")) else NULL,
    composition = if (repetition == 1L) composition else NULL,
    truth = if (repetition == 1L) data.frame(z = grid, value = truth_centered) else NULL,
    sample_focal = if (repetition == 1L) data$x else NULL
  )
}

run_s23 <- function(repetitions = 300L, n_evaluation = 1000L,
                    sigma = 0.45, alpha = 0.05, epsilon = 0.10,
                    grid = seq(-0.75, 0.75, length.out = 21L),
                    beta = 1.5, strength = 4, seed = 20260925L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  grid <- assert_grid(grid)
  evaluate <- function(index) s23_one(
    index, n_evaluation, sigma, alpha, epsilon, grid, beta, strength, seed
  )
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(repetitions), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(repetitions), evaluate)
  }
  list(
    metrics = do.call(rbind, lapply(pieces, `[[`, "metrics")),
    curves = pieces[[1]]$curves,
    composition = pieces[[1]]$composition,
    truth = pieces[[1]]$truth,
    sample_focal = pieces[[1]]$sample_focal,
    settings = list(repetitions = repetitions, n_evaluation = n_evaluation,
                    sigma = sigma, alpha = alpha, epsilon = epsilon,
                    beta = beta, strength = strength, grid = grid)
  )
}

summarize_s23 <- function(simulation) {
  groups <- split(
    simulation$metrics,
    interaction(simulation$metrics$extension, simulation$metrics$method,
                drop = TRUE)
  )
  do.call(rbind, lapply(groups, function(part) {
    data.frame(
      extension = part$extension[[1]], method = part$method[[1]],
      normalized_shape_error = mean(part$normalized_shape_error),
      normalized_shape_error_mc_se = stats::sd(part$normalized_shape_error) /
        sqrt(nrow(part)),
      rmse = mean(part$rmse),
      extension_sensitivity = mean(part$extension_sensitivity),
      extension_sensitivity_mc_se = stats::sd(part$extension_sensitivity) /
        sqrt(nrow(part)),
      gamma = mean(part$gamma), mean_isc = mean(part$mean_isc),
      rcpd_min_coverage = mean(part$rcpd_min_coverage),
      rcpd_simultaneous_mass = mean(part$rcpd_simultaneous_mass),
      rcpd_ess_fraction = mean(part$rcpd_ess_fraction),
      stringsAsFactors = FALSE
    )
  }))
}

plot_s23 <- function(simulation,
                     path = file.path("results", "s23-composition-trap.png")) {
  method_order <- c("Truth", "PDP", "Pointwise trim", "Hard CSPD", "RCPD", "ALE")
  colors <- c(
    Truth = pdp_palette[["ink"]], PDP = pdp_palette[["orange"]],
    `Pointwise trim` = pdp_palette[["purple"]],
    `Hard CSPD` = pdp_palette[["green"]], RCPD = pdp_palette[["blue"]],
    ALE = pdp_palette[["gold"]]
  )
  curves <- simulation$curves[
    simulation$curves$extension == "Upward extension", , drop = FALSE
  ]
  truth <- transform(simulation$truth, method = "Truth",
                     extension = "Upward extension")
  curves <- rbind(truth, curves)
  curves$method <- factor(curves$method, levels = method_order)
  rug <- data.frame(z = simulation$sample_focal)
  rug <- rug[rug$z >= min(curves$z) & rug$z <= max(curves$z), , drop = FALSE]
  p_curves <- ggplot2::ggplot(
    curves, ggplot2::aes(z, value, color = method, linetype = method)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_rug(
      data = rug, ggplot2::aes(x = z), inherit.aes = FALSE, sides = "b",
      alpha = 0.13, color = pdp_palette[["mid_grey"]]
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash",
                                              "longdash", "twodash", "dotted")) +
    ggplot2::labs(x = expression(z), y = "Centered feature effect",
                  color = NULL, linetype = NULL) +
    theme_pdp()

  composition <- rbind(
    data.frame(z = simulation$composition$z,
               value = simulation$composition$pointwise_mean_u,
               reference = "Pointwise-trimmed"),
    data.frame(z = simulation$composition$z,
               value = simulation$composition$original_mean_u,
               reference = "Original P"),
    data.frame(z = simulation$composition$z,
               value = simulation$composition$hard_mean_u,
               reference = "Hard common"),
    data.frame(z = simulation$composition$z,
               value = simulation$composition$relaxed_mean_u,
               reference = "RCPD")
  )
  p_composition <- ggplot2::ggplot(
    composition, ggplot2::aes(z, value, color = reference, linetype = reference)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = pdp_palette[["light_grey"]]) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::scale_color_manual(values = c(
      `Pointwise-trimmed` = pdp_palette[["purple"]],
      `Original P` = pdp_palette[["orange"]],
      `Hard common` = pdp_palette[["green"]], RCPD = pdp_palette[["blue"]]
    )) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
    ggplot2::labs(x = expression(z), y = "Mean context U",
                  color = NULL, linetype = NULL) +
    theme_pdp()

  summary <- summarize_s23(simulation)
  upward <- summary[summary$extension == "Upward extension", ]
  upward$method <- factor(upward$method, levels = method_order[-1])
  p_error <- ggplot2::ggplot(
    upward, ggplot2::aes(method, normalized_shape_error, color = method)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, normalized_shape_error - 2 * normalized_shape_error_mc_se),
                   ymax = normalized_shape_error + 2 * normalized_shape_error_mc_se),
      width = 0.15, linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "Normalized shape error") +
    theme_pdp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "none")

  sensitivity <- summary[!duplicated(summary$method), ]
  sensitivity$method <- factor(sensitivity$method, levels = method_order[-1])
  p_sensitivity <- ggplot2::ggplot(
    sensitivity, ggplot2::aes(method, extension_sensitivity, color = method)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, extension_sensitivity - 2 * extension_sensitivity_mc_se),
                   ymax = extension_sensitivity + 2 * extension_sensitivity_mc_se),
      width = 0.15, linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = NULL, y = "Sensitivity to off-support extension") +
    theme_pdp() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "none")

  save_pdp_panel_set(
    list(p_curves, p_composition, p_error, p_sensitivity), path,
    widths = 5.2, heights = 4.0
  )
}
