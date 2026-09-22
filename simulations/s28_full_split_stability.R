# Full-pipeline repeated-split stability for the held-out applications.

s28_extract <- function(case, application = NULL) {
  if (is.null(application)) application <- unique(case$dataset)
  model_groups <- if ("model_case" %in% names(case)) unique(case$model_case) else
    unique(case$model)
  do.call(rbind, lapply(model_groups, function(model_label) {
    part <- if ("model_case" %in% names(case))
      case[case$model_case == model_label, , drop = FALSE] else case
    pd_centered <- part$pdp - mean(part$pdp)
    hard_centered <- part$cspd - mean(part$cspd)
    relaxed_centered <- part$relaxed_cspd - mean(part$relaxed_cspd)
    difference <- relaxed_centered - pd_centered
    maximum <- which.max(abs(difference))
    amplitude <- diff(range(pd_centered))
    data.frame(
      application = application, model = model_label,
      pd_amplitude = amplitude,
      hard_distance = max(abs(hard_centered - pd_centered)),
      relaxed_distance = max(abs(difference)),
      normalized_hard_distance = if (amplitude > 0)
        max(abs(hard_centered - pd_centered)) / amplitude else NA_real_,
      normalized_relaxed_distance = if (amplitude > 0)
        max(abs(difference)) / amplitude else NA_real_,
      max_difference_z = part$z[[maximum]],
      max_difference_sign = sign(difference[[maximum]]),
      gamma = unique(part$gamma), n_common = unique(part$n_common),
      ess = unique(part$relaxed_ess), kl = unique(part$relaxed_kl),
      min_coverage = unique(part$relaxed_min_coverage),
      simultaneous_valid_mass = unique(part$relaxed_simultaneous_valid_mass),
      predictive_loss = if ("heldout_mae" %in% names(part))
        unique(part$heldout_mae) else unique(part$heldout_brier),
      stringsAsFactors = FALSE
    )
  }))
}

run_s28 <- function(wine_repetitions = 50L, bank_repetitions = 20L,
                    shopper_repetitions = 20L, seed = 20260930L,
                    workers = if (.Platform$OS.type == "windows") 1L else 4L) {
  design <- rbind(
    data.frame(application = "Wine", repetition = seq_len(wine_repetitions)),
    data.frame(application = "Bank marketing", repetition = seq_len(bank_repetitions)),
    data.frame(application = "Online shoppers", repetition = seq_len(shopper_repetitions))
  )
  evaluate <- function(index) {
    setting <- design[index, ]
    case_seed <- seed + index
    case <- switch(
      setting$application,
      Wine = run_s16(seed = case_seed),
      `Bank marketing` = run_s17_bank(seed = case_seed),
      `Online shoppers` = run_s17_shoppers(seed = case_seed)
    )
    output <- s28_extract(case, setting$application)
    output$repetition <- setting$repetition
    output
  }
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(design)), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(nrow(design)), evaluate)
  }
  do.call(rbind, pieces)
}

summarize_s28 <- function(simulation) {
  key <- interaction(simulation$application, simulation$model, drop = TRUE)
  pieces <- split(simulation, key)
  do.call(rbind, lapply(pieces, function(part) {
    q <- function(x, p) stats::quantile(x, p, na.rm = TRUE, names = FALSE)
    positive <- mean(part$max_difference_sign > 0, na.rm = TRUE)
    data.frame(
      application = part$application[[1]], model = part$model[[1]],
      repetitions = nrow(part),
      median_normalized_distance = stats::median(
        part$normalized_relaxed_distance, na.rm = TRUE
      ),
      distance_05 = q(part$normalized_relaxed_distance, 0.05),
      distance_95 = q(part$normalized_relaxed_distance, 0.95),
      dominant_sign_rate = max(positive, 1 - positive),
      median_gamma = stats::median(part$gamma),
      gamma_05 = q(part$gamma, 0.05), gamma_95 = q(part$gamma, 0.95),
      median_ess = stats::median(part$ess),
      median_simultaneous_mass = stats::median(part$simultaneous_valid_mass),
      median_predictive_loss = stats::median(part$predictive_loss),
      stringsAsFactors = FALSE
    )
  }))
}

plot_s28 <- function(simulation,
                     path = file.path("results", "s28-full-split-stability.png")) {
  simulation$case <- paste(simulation$application, simulation$model, sep = ": ")
  p_distance <- ggplot2::ggplot(
    simulation, ggplot2::aes(normalized_relaxed_distance, case,
                             color = application)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25, width = 0.58,
                          linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c(
      `Wine` = pdp_palette[["purple"]],
      `Bank marketing` = pdp_palette[["blue"]],
      `Online shoppers` = pdp_palette[["orange"]]
    )) +
    ggplot2::labs(x = "Normalized RCPD-PDP shape distance", y = NULL,
                  color = NULL) +
    theme_pdp() + ggplot2::theme(legend.position = "none")

  p_gamma <- ggplot2::ggplot(
    simulation, ggplot2::aes(gamma, case, color = application)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25, width = 0.58,
                          linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c(
      `Wine` = pdp_palette[["purple"]],
      `Bank marketing` = pdp_palette[["blue"]],
      `Online shoppers` = pdp_palette[["orange"]]
    )) +
    ggplot2::labs(x = "Hard common mass", y = NULL, color = NULL) +
    theme_pdp() + ggplot2::theme(legend.position = "none")
  save_pdp_panel_set(list(p_distance, p_gamma), path,
                     widths = 5.2, heights = 4.0)
}
