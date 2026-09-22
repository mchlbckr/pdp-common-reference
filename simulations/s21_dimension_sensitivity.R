# Simulation S21: context dimension and conditional-versus-joint validity.
# An equicorrelated Gaussian design represents a one-factor structure with the
# same pairwise correlation at every dimension.

simulate_s21_data <- function(n, context_dimension, rho, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dimension <- context_dimension + 1L
  covariance <- matrix(rho, nrow = dimension, ncol = dimension)
  diag(covariance) <- 1
  values <- matrix(stats::rnorm(n * dimension), nrow = n) %*% chol(covariance)
  output <- as.data.frame(values)
  names(output) <- c("x1", paste0("u", seq_len(context_dimension)))
  attr(output, "covariance") <- covariance
  output
}

s21_oracle_validity <- function(context_dimension, rho, alpha, score) {
  dimension <- context_dimension + 1L
  covariance <- matrix(rho, nrow = dimension, ncol = dimension)
  diag(covariance) <- 1
  if (score == "conditional") {
    context_covariance <- covariance[-1, -1, drop = FALSE]
    beta <- as.numeric(covariance[1, -1, drop = FALSE] %*% solve(context_covariance))
    residual_sd <- as.numeric(sqrt(covariance[1, 1] -
                                     covariance[1, -1, drop = FALSE] %*%
                                       solve(context_covariance, covariance[-1, 1, drop = FALSE])))
    cutoff <- stats::qnorm(1 - alpha / 2)
    return(function(newdata) {
      context <- as.matrix(newdata[paste0("u", seq_len(context_dimension))])
      abs((newdata$x1 - as.numeric(context %*% beta)) / residual_sd) <= cutoff
    })
  }
  if (score == "joint") {
    precision <- solve(covariance)
    cutoff <- stats::qchisq(1 - alpha, df = dimension)
    return(function(newdata) {
      values <- as.matrix(newdata)
      rowSums((values %*% precision) * values) <= cutoff
    })
  }
  stop("score must be 'conditional' or 'joint'.", call. = FALSE)
}

run_s21 <- function(context_dimensions = c(1L, 4L, 9L), repetitions = 100L,
                    n_evaluation = 1000L, rho = 0.5, alpha = 0.05,
                    epsilon = 0.05, grid = seq(-1.1, 1.1, length.out = 21L),
                    seed = 20260924L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  design <- expand.grid(
    context_dimension = context_dimensions,
    score = c("conditional", "joint"), repetition = seq_len(repetitions),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  evaluate <- function(index) {
    setting <- design[index, ]
    data <- simulate_s21_data(
      n_evaluation, setting$context_dimension, rho, seed = seed + index
    )
    valid <- s21_oracle_validity(
      setting$context_dimension, rho, alpha, setting$score
    )
    validity <- validity_matrix(data, "x1", grid, valid)
    diagnostics <- validity_pattern_diagnostics(validity)
    fitted_epsilon <- max(epsilon, diagnostics$epsilon_min + 1e-6)
    reference <- fit_relaxed_reference(validity, epsilon = fitted_epsilon, tolerance = 2e-5)
    data.frame(
      context_dimension = setting$context_dimension,
      score = setting$score, repetition = setting$repetition,
      mean_isc = mean(validity), gamma = mean(rowSums(validity) == ncol(validity)),
      epsilon_min = diagnostics$epsilon_min,
      fitted_epsilon = fitted_epsilon,
      effective_sample_fraction = reference$effective_sample_size / n_evaluation,
      observed_patterns = diagnostics$observed_patterns,
      singleton_pattern_fraction = diagnostics$singleton_pattern_fraction,
      stringsAsFactors = FALSE
    )
  }
  indices <- seq_len(nrow(design))
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(indices, evaluate, mc.cores = workers)
  } else {
    lapply(indices, evaluate)
  }
  do.call(rbind, pieces)
}

summarize_s21 <- function(simulation) {
  measures <- c("mean_isc", "gamma", "epsilon_min", "fitted_epsilon",
                "effective_sample_fraction", "observed_patterns",
                "singleton_pattern_fraction")
  pieces <- split(simulation, interaction(simulation$context_dimension,
                                          simulation$score, drop = TRUE))
  do.call(rbind, lapply(pieces, function(part) {
    output <- part[1, c("context_dimension", "score"), drop = FALSE]
    for (measure in measures) {
      output[[measure]] <- mean(part[[measure]])
      output[[paste0(measure, "_mc_se")]] <- stats::sd(part[[measure]]) / sqrt(nrow(part))
    }
    output
  }))
}

plot_s21 <- function(summary, path = file.path("results", "s21-dimension-sensitivity.png")) {
  labels <- c(conditional = "Conditional score", joint = "Joint score")
  summary$display_score <- factor(labels[summary$score], levels = labels)
  colors <- c("Conditional score" = pdp_palette[["blue"]],
              "Joint score" = pdp_palette[["orange"]])
  make_panel <- function(measure, y_label, limits = NULL) {
    plot <- ggplot2::ggplot(
      summary, ggplot2::aes(context_dimension, .data[[measure]],
                            color = display_score, linetype = display_score)
    ) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = sort(unique(summary$context_dimension))) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
      ggplot2::labs(x = "Context dimension", y = y_label, color = NULL, linetype = NULL) +
      theme_pdp()
    if (!is.null(limits)) {
      plot <- plot + ggplot2::coord_cartesian(ylim = limits)
    }
    plot
  }
  epsilon_panel <- make_panel(
    "epsilon_min", expression(hat(epsilon)[min]), c(-0.005, 0.05)
  ) +
    ggplot2::annotate(
      "text", x = max(summary$context_dimension), y = 0.035,
      label = "All values = 0", hjust = 1, size = 3,
      color = pdp_palette[["mid_grey"]]
    )
  panels <- list(
    make_panel("mean_isc", expression(widehat(ISC)^P(G)), c(0, 1)),
    make_panel("gamma", expression(hat(gamma)(G)), c(0, 1)),
    epsilon_panel,
    make_panel("effective_sample_fraction", "Effective sample fraction", c(0, 1))
  )
  save_pdp_panel_set(panels, path, widths = 5.2, heights = 4.0)
}
