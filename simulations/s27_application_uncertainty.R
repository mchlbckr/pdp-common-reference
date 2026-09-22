# Evaluation-resampling stability bands for the fixed fitted application models.

run_s27 <- function(wine, large_cases, repetitions = 200L, level = 0.95,
                    seed = 20260929L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  fits <- list()
  wine_validity <- attr(wine, "evaluation_validity")
  wine_predictions <- attr(wine, "prediction_matrices")
  for (model_name in names(wine_predictions)) {
    label <- unique(wine$model_case[wine$model_name == model_name])
    fits[[paste("Wine", label, sep = ": ")]] <-
      bootstrap_reference_curve_stability(
        wine_validity, wine_predictions[[model_name]],
        grid = unique(wine$z), epsilon = unique(wine$relaxed_epsilon),
        repetitions = repetitions, level = level,
        seed = seed + length(fits), workers = workers
      )
  }
  for (case in large_cases) {
    label <- unique(case$dataset)
    fits[[label]] <- bootstrap_reference_curve_stability(
      attr(case, "evaluation_validity"), attr(case, "prediction_matrix"),
      grid = case$z, epsilon = unique(case$relaxed_epsilon),
      repetitions = repetitions, level = level,
      seed = seed + length(fits), workers = workers
    )
  }
  differences <- do.call(rbind, lapply(names(fits), function(application) {
    output <- fits[[application]]$differences
    output$application <- application
    output
  }))
  bands <- do.call(rbind, lapply(names(fits), function(application) {
    output <- fits[[application]]$bands
    output$application <- application
    output
  }))
  metrics <- do.call(rbind, lapply(names(fits), function(application) {
    output <- fits[[application]]$metrics
    output$application <- application
    output
  }))
  list(fits = fits, differences = differences, bands = bands,
       metrics = metrics, repetitions = repetitions, level = level)
}

summarize_s27 <- function(simulation) {
  pieces <- split(simulation$metrics, simulation$metrics$application)
  do.call(rbind, lapply(pieces, function(part) {
    quantile_value <- function(x, probability) stats::quantile(
      x, probability, na.rm = TRUE, names = FALSE
    )
    data.frame(
      application = part$application[[1]],
      median_gamma = stats::median(part$gamma),
      gamma_05 = quantile_value(part$gamma, 0.05),
      gamma_95 = quantile_value(part$gamma, 0.95),
      median_ess = stats::median(part$ess),
      median_simultaneous_mass = stats::median(part$simultaneous_valid_mass),
      median_normalized_hard_distance = stats::median(
        part$normalized_hard_distance, na.rm = TRUE
      ),
      median_normalized_relaxed_distance = stats::median(
        part$normalized_relaxed_distance, na.rm = TRUE
      ),
      relaxed_distance_05 = quantile_value(
        part$normalized_relaxed_distance, 0.05
      ),
      relaxed_distance_95 = quantile_value(
        part$normalized_relaxed_distance, 0.95
      ),
      stringsAsFactors = FALSE
    )
  }))
}

plot_s27 <- function(simulation,
                     path = file.path("results", "s27-application-uncertainty.png")) {
  differences <- simulation$differences[
    simulation$differences$contrast == "RCPD minus PDP", , drop = FALSE
  ]
  p_difference <- ggplot2::ggplot(
    differences, ggplot2::aes(z, estimate)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = pdp_palette[["mid_grey"]],
                        linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = pdp_palette[["sky"]], alpha = 0.24
    ) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.75) +
    ggplot2::facet_wrap(~application, scales = "free", ncol = 2) +
    ggplot2::labs(x = "Focal-feature grid", y = "Centered RCPD minus PDP") +
    theme_pdp()

  metrics <- simulation$metrics
  p_distance <- ggplot2::ggplot(
    metrics, ggplot2::aes(normalized_relaxed_distance, application,
                          color = application)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c(
      pdp_palette[["blue"]], pdp_palette[["orange"]],
      pdp_palette[["green"]], pdp_palette[["purple"]]
    )) +
    ggplot2::labs(x = "Normalized RCPD-PDP shape distance",
                  y = NULL) +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  save_pdp_panel_set(
    list(p_difference, p_distance), path,
    widths = 10, heights = c(5.5, 3.8)
  )
}
