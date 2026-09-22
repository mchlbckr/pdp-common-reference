# Simulation S11: oracle-target bootstrap calibration across sample sizes.

run_s11 <- function(
  split_sizes = c(75L, 150L, 300L),
  repetitions = 60L,
  bootstrap_repetitions = 50L,
  rho = 0.7,
  alpha = 0.05,
  grid = seq(-0.75, 0.75, length.out = 11L),
  seed = 20260919L
) {
  studies <- lapply(seq_along(split_sizes), function(index) {
    split_size <- split_sizes[[index]]
    simulation <- run_s10(
      n_training = split_size,
      n_calibration = split_size,
      n_evaluation = 4L * split_size,
      repetitions = repetitions,
      bootstrap_repetitions = bootstrap_repetitions,
      rho = rho, alpha = alpha, grid = grid, seed = seed + index
    )
    data.frame(
      n_training = split_size,
      n_calibration = split_size,
      n_evaluation = 4L * split_size,
      total_n = 6L * split_size,
      wald_coverage = mean(simulation$covered_wald),
      bootstrap_coverage = mean(simulation$covered_bootstrap),
      bootstrap_width = mean(simulation$bootstrap_width)
    )
  })
  do.call(rbind, studies)
}

plot_s11 <- function(simulation, path = file.path("results", "s11-bootstrap-sample-size.png")) {
  coverage <- rbind(
    data.frame(total_n = simulation$total_n, coverage = simulation$wald_coverage,
               method = "Plug-in Wald"),
    data.frame(total_n = simulation$total_n, coverage = simulation$bootstrap_coverage,
               method = "Outer bootstrap")
  )
  coverage$method <- factor(coverage$method,
                            levels = c("Plug-in Wald", "Outer bootstrap"))
  colors <- c("Plug-in Wald" = pdp_palette[["orange"]],
              "Outer bootstrap" = pdp_palette[["blue"]])
  p_coverage <- ggplot2::ggplot(
    coverage, ggplot2::aes(total_n, coverage, color = method,
                           shape = method, linetype = method)
  ) +
    ggplot2::geom_hline(yintercept = 0.95, color = pdp_palette[["ink"]],
                        linetype = "dotted", linewidth = 0.55) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.2, stroke = 0.7) +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = "Nominal 95%",
                      hjust = 1.05, vjust = -0.5, size = 3) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_shape_manual(values = c("Plug-in Wald" = 17,
                                           "Outer bootstrap" = 16)) +
    ggplot2::scale_linetype_manual(values = c("Plug-in Wald" = "dashed",
                                              "Outer bootstrap" = "solid")) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "Total split-sample size", y = "Mean oracle-target coverage",
                  color = NULL, shape = NULL, linetype = NULL) +
    theme_pdp()
  p_width <- ggplot2::ggplot(simulation, ggplot2::aes(total_n, bootstrap_width)) +
    ggplot2::geom_line(color = pdp_palette[["blue"]], linewidth = 0.7) +
    ggplot2::geom_point(color = pdp_palette[["blue"]], size = 2.2) +
    ggplot2::labs(x = "Total split-sample size", y = "Mean percentile-interval width") +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  save_pdp_panel_set(list(p_coverage, p_width), path,
                     widths = 5.2, heights = 4.0)
}
