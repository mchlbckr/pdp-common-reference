# Simulation S26: finite display grids versus the exact empirical
# semi-infinite interval-validity projection.

s26_one <- function(repetition, n_evaluation, rho, alpha, epsilon, interval,
                    grid_sizes, validation_size, seed) {
  data <- simulate_s1_data(n_evaluation, rho, seed + repetition)
  valid <- s14_oracle_validity(rho, alpha)
  cutoff <- stats::qnorm(1 - alpha / 2)
  residual_scale <- sqrt(1 - rho^2)
  lower <- rho * data$x2 - cutoff * residual_scale
  upper <- rho * data$x2 + cutoff * residual_scale
  validation_grid <- seq(interval[[1]], interval[[2]],
                         length.out = validation_size)
  exact <- fit_adaptive_interval_relaxed_reference(
    lower, upper, interval, epsilon = epsilon,
    initial_grid = c(interval[[1]], mean(interval), interval[[2]]),
    tolerance = 2e-5
  )
  if (!exact$converged) {
    stop("exact interval exchange solver did not converge.", call. = FALSE)
  }
  exact_curve <- estimate_weighted_pdp(
    data, "x1", validation_grid, s18_prediction, exact$weights
  )$weighted_pdp
  critical_h <- interval_validity_matrix(lower, upper, exact$critical_grid)

  fixed <- do.call(rbind, lapply(grid_sizes, function(k) {
    grid <- seq(interval[[1]], interval[[2]], length.out = k)
    h <- validity_matrix(data, "x1", grid, valid)
    fit <- fit_relaxed_reference(h, epsilon = epsilon, tolerance = 2e-5)
    continuum_coverage <- as.numeric(crossprod(
      fit$weights, unclass(critical_h) * 1
    ))
    curve <- estimate_weighted_pdp(
      data, "x1", validation_grid, s18_prediction, fit$weights
    )$weighted_pdp
    data.frame(
      repetition = repetition, method = "Fixed grid", grid_size = k,
      hard_gamma = mean(rowSums(h) == ncol(h)),
      min_continuum_coverage = min(continuum_coverage),
      mean_continuum_coverage = mean(continuum_coverage),
      curve_distance_to_exact = centered_curve_distance(curve, exact_curve),
      kl = fit$kl, ess_fraction = fit$effective_sample_size / n_evaluation,
      constraints_used = k, active_constraints =
        fit$pattern_diagnostics$active_constraints,
      stringsAsFactors = FALSE
    )
  }))

  exact_row <- data.frame(
    repetition = repetition, method = "Exact exchange",
    grid_size = validation_size,
    hard_gamma = mean(rowSums(critical_h) == ncol(critical_h)),
    min_continuum_coverage = exact$min_critical_coverage,
    mean_continuum_coverage = exact$mean_critical_coverage,
    curve_distance_to_exact = 0,
    kl = exact$kl,
    ess_fraction = exact$effective_sample_size / n_evaluation,
    constraints_used = length(exact$active_grid),
    active_constraints = exact$fit$pattern_diagnostics$active_constraints,
    stringsAsFactors = FALSE
  )
  rbind(fixed, exact_row)
}

run_s26 <- function(repetitions = 100L, n_evaluation = 1000L,
                    rho = 0.8, alpha = 0.05, epsilon = 0.10,
                    interval = c(-1.1, 1.1),
                    grid_sizes = c(3L, 5L, 11L, 21L, 41L, 81L),
                    validation_size = 161L, seed = 20260928L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  if (any(grid_sizes < 2L) || any(grid_sizes > validation_size)) {
    stop("grid_sizes must lie between 2 and validation_size.", call. = FALSE)
  }
  evaluate <- function(index) s26_one(
    index, n_evaluation, rho, alpha, epsilon, interval, grid_sizes,
    validation_size, seed
  )
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(repetitions), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(repetitions), evaluate)
  }
  output <- do.call(rbind, pieces)
  attr(output, "settings") <- list(
    repetitions = repetitions, n_evaluation = n_evaluation, rho = rho,
    alpha = alpha, epsilon = epsilon, interval = interval,
    grid_sizes = grid_sizes, validation_size = validation_size
  )
  output
}

summarize_s26 <- function(simulation) {
  key <- interaction(simulation$method, simulation$grid_size, drop = TRUE)
  pieces <- split(simulation, key)
  do.call(rbind, lapply(pieces, function(part) {
    measures <- c("hard_gamma", "min_continuum_coverage",
                  "mean_continuum_coverage", "curve_distance_to_exact",
                  "kl", "ess_fraction", "constraints_used",
                  "active_constraints")
    output <- data.frame(method = part$method[[1]],
                         grid_size = part$grid_size[[1]])
    for (measure in measures) {
      output[[measure]] <- mean(part[[measure]])
      output[[paste0(measure, "_mc_se")]] <- stats::sd(part[[measure]]) /
        sqrt(nrow(part))
    }
    output
  }))
}

plot_s26 <- function(summary,
                     path = file.path("results", "s26-grid-resolution.png")) {
  fixed <- summary[summary$method == "Fixed grid", ]
  adaptive <- summary[summary$method == "Exact exchange", ]
  epsilon <- 0.10
  p_gamma <- ggplot2::ggplot(fixed, ggplot2::aes(grid_size, hard_gamma)) +
    ggplot2::geom_line(color = pdp_palette[["green"]], linewidth = 0.75) +
    ggplot2::geom_point(color = pdp_palette[["green"]], size = 2.1) +
    ggplot2::scale_x_log10(breaks = fixed$grid_size) +
    ggplot2::labs(x = "Grid points K", y = "Hard common mass") +
    theme_pdp()

  p_coverage <- ggplot2::ggplot(
    fixed, ggplot2::aes(grid_size, min_continuum_coverage)
  ) +
    ggplot2::geom_hline(yintercept = 1 - epsilon, linetype = "dashed",
                        color = pdp_palette[["ink"]]) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.75) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.1) +
    ggplot2::geom_point(
      data = adaptive,
      ggplot2::aes(x = max(fixed$grid_size), y = min_continuum_coverage),
      inherit.aes = FALSE, shape = 17, size = 2.5,
      color = pdp_palette[["orange"]]
    ) +
    ggplot2::scale_x_log10(breaks = fixed$grid_size) +
    ggplot2::labs(x = "Display-grid points K", y = "Minimum continuum coverage") +
    theme_pdp()

  p_distance <- ggplot2::ggplot(
    fixed, ggplot2::aes(grid_size, curve_distance_to_exact)
  ) +
    ggplot2::geom_line(color = pdp_palette[["purple"]], linewidth = 0.75) +
    ggplot2::geom_point(color = pdp_palette[["purple"]], size = 2.1) +
    ggplot2::geom_point(
      data = adaptive,
      ggplot2::aes(x = max(fixed$grid_size), y = curve_distance_to_exact),
      inherit.aes = FALSE, shape = 17, size = 2.5,
      color = pdp_palette[["orange"]]
    ) +
    ggplot2::scale_x_log10(breaks = fixed$grid_size) +
    ggplot2::labs(x = "Display-grid points K", y = "Shape distance to exact solution") +
    theme_pdp()

  constraints <- rbind(
    data.frame(grid_size = fixed$grid_size,
               constraints = fixed$constraints_used,
               method = "Fixed grid"),
    data.frame(grid_size = max(fixed$grid_size),
               constraints = adaptive$constraints_used,
               method = "Exact exchange")
  )
  p_constraints <- ggplot2::ggplot(
    constraints, ggplot2::aes(grid_size, constraints, color = method,
                               shape = method, linetype = method)
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.1) +
    ggplot2::scale_x_log10(breaks = fixed$grid_size) +
    ggplot2::scale_color_manual(values = c(
      `Fixed grid` = pdp_palette[["mid_grey"]],
      `Exact exchange` = pdp_palette[["orange"]]
    )) +
    ggplot2::scale_shape_manual(values = c(`Fixed grid` = 16,
                                            `Exact exchange` = 17)) +
    ggplot2::scale_linetype_manual(values = c(`Fixed grid` = "solid",
                                               `Exact exchange` = "blank")) +
    ggplot2::labs(x = "Grid points K", y = "Constraints used",
                  color = NULL, shape = NULL, linetype = NULL) +
    theme_pdp()

  save_pdp_panel_set(
    list(p_gamma, p_coverage, p_distance, p_constraints), path,
    widths = 5.2, heights = 4.0
  )
}
