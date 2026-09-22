# Shared publication theme and color system for the main-paper figures.

pdp_palette <- c(
  ink = "#222222",
  blue = "#0072B2",
  orange = "#D55E00",
  green = "#009E73",
  purple = "#CC79A7",
  sky = "#56B4E9",
  gold = "#E69F00",
  light_grey = "#D9D9D9",
  mid_grey = "#7A7A7A"
)

theme_pdp <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.3),
      axis.line = ggplot2::element_line(color = pdp_palette[["ink"]], linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = pdp_palette[["ink"]], linewidth = 0.35),
      axis.title = ggplot2::element_text(color = pdp_palette[["ink"]]),
      axis.text = ggplot2::element_text(color = pdp_palette[["ink"]]),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      plot.subtitle = ggplot2::element_text(color = pdp_palette[["mid_grey"]],
                                            size = ggplot2::rel(0.9)),
      plot.tag = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(5, 7, 5, 5)
    )
}

save_pdp_plot <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, bg = "white")
  path
}

# Save plots as separate files so that Quarto, rather than an image compositor,
# owns subfigure labels, captions, layout, and panel-level cross-references.
save_pdp_panel_set <- function(plots, path, widths = 5.2, heights = 4.0,
                               dpi = 300) {
  if (!is.list(plots) || length(plots) < 2L) {
    stop("plots must be a list containing at least two panels.", call. = FALSE)
  }
  extension <- tools::file_ext(path)
  if (!nzchar(extension)) extension <- "png"
  stem <- sub(paste0("\\.", extension, "$"), "", basename(path))
  paths <- file.path(
    dirname(path), paste0(stem, "-", letters[seq_along(plots)], ".", extension)
  )
  widths <- rep(widths, length.out = length(plots))
  heights <- rep(heights, length.out = length(plots))
  for (index in seq_along(plots)) {
    save_pdp_plot(
      plots[[index]], paths[[index]], widths[[index]], heights[[index]], dpi
    )
  }
  paths
}
