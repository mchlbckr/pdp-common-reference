# Simulation S14: common protocol across models, grids, learners, and baselines.

s14_curve_amplitude <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else diff(range(x))
}

s14_shape_distance <- function(x, target) {
  if (length(x) != length(target) || any(!is.finite(x)) || any(!is.finite(target))) NA_real_
  else centered_curve_distance(x, target)
}

s14_fit_model <- function(training_data, model_class) {
  if (model_class == "linear") {
    return(stats::lm(y ~ x1 * x2 + I(x2^2), data = training_data))
  }
  if (model_class == "tree") {
    return(rpart::rpart(y ~ x1 + x2, data = training_data,
                        control = rpart::rpart.control(cp = 0.01, minsplit = 20L)))
  }
  if (model_class == "gam") {
    return(mgcv::bam(y ~ s(x1, k = 6) + s(x2, k = 6) + ti(x1, x2, k = c(4, 4)),
                     data = training_data, method = "fREML", discrete = TRUE))
  }
  if (model_class == "neural_net") {
    return(nnet::nnet(y ~ x1 + x2, data = training_data, size = 6, linout = TRUE,
                      decay = 0.01, maxit = 300, trace = FALSE))
  }
  if (model_class == "random_forest") {
    return(randomForest::randomForest(y ~ x1 + x2, data = training_data,
                                      ntree = 150L, mtry = 2L, nodesize = 5L))
  }
  if (model_class == "boosted_tree") {
    return(gbm::gbm(y ~ x1 + x2, data = training_data, distribution = "gaussian",
                    n.trees = 150L, interaction.depth = 2L, shrinkage = 0.05,
                    n.minobsinnode = 10L, verbose = FALSE))
  }
  stop("unknown model_class", call. = FALSE)
}

s14_predict <- function(model) {
  if (inherits(model, "gbm")) {
    return(function(newdata) as.numeric(stats::predict(model, newdata = newdata, n.trees = model$n.trees)))
  }
  function(newdata) as.numeric(stats::predict(model, newdata = newdata))
}

s14_add_response <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data$y <- sin(data$x1) + 0.5 * data$x2^2 + 0.9 * data$x1 * data$x2 + stats::rnorm(nrow(data), sd = 0.15)
  data
}

s14_oracle_validity <- function(rho, alpha) {
  cutoff <- stats::qnorm(1 - alpha / 2)
  function(newdata) {
    residual <- (newdata$x1 - rho * newdata$x2) / sqrt(1 - rho^2)
    abs(residual) <= cutoff
  }
}

s14_one <- function(rho, grid_half_width, model_class, repetition,
                    n_training, n_calibration, n_evaluation, alpha, k, seed) {
  set.seed(seed)
  training <- s14_add_response(simulate_s1_data(n_training, rho))
  calibration <- simulate_s1_data(n_calibration, rho)
  evaluation <- simulate_s1_data(n_evaluation, rho)
  model <- s14_fit_model(training, model_class)
  prediction <- s14_predict(model)
  grid <- seq(-grid_half_width, grid_half_width, length.out = 21L)
  oracle <- s14_oracle_validity(rho, alpha)
  conformal_fit <- fit_split_conformal_validity(training[, c("x1", "x2")], calibration, "x1", alpha)
  density_fit <- fit_conditional_density_validity(training[, c("x1", "x2")], calibration, "x1", alpha)
  local_scale_fit <- fit_locally_scaled_validity(training[, c("x1", "x2")], calibration, "x1", alpha)
  knn_fit <- fit_knn_validity(training[, c("x1", "x2")], k = k, alpha = alpha,
                              calibration_data = calibration)
  learners <- list(
    oracle = oracle,
    split_conformal = function(newdata) predict_validity(conformal_fit, newdata),
    conditional_density = function(newdata) predict_conditional_density_validity(density_fit, newdata),
    locally_scaled = function(newdata) predict_locally_scaled_validity(local_scale_fit, newdata),
    knn = function(newdata) predict_knn_validity(knn_fit, newdata)
  )
  calibration_exceedance <- c(
    oracle = NA_real_, split_conformal = conformal_fit$calibration_exceedance,
    conditional_density = density_fit$calibration_exceedance,
    locally_scaled = local_scale_fit$calibration_exceedance,
    knn = knn_fit$calibration_exceedance
  )
  residual_validity <- validity_matrix(evaluation, "x1", grid, learners$split_conformal)
  density_validity <- validity_matrix(evaluation, "x1", grid, learners$conditional_density)
  residual_density_identical <- identical(residual_validity, density_validity)
  if (!residual_density_identical) {
    stop("split-residual and homoskedastic conditional-density validity matrices differ.",
         call. = FALSE)
  }

  # Baselines not indexed by a validity learner: their estimands are fixed by
  # their own reference rules.
  pdp <- estimate_pdp(evaluation, "x1", grid, prediction)
  conditional <- estimate_conditional_curve(evaluation, "x1", grid, prediction, bandwidth = 0.2)
  breaks <- unique(stats::quantile(evaluation$x1, probs = seq(0, 1, length.out = 13L)))
  ale <- estimate_ale(evaluation, "x1", breaks, prediction)
  low_subgroup <- estimate_subgroup_pdp(evaluation, "x1", grid, prediction,
                                        function(data) data$x2 <= stats::median(evaluation$x2))
  high_subgroup <- estimate_subgroup_pdp(evaluation, "x1", grid, prediction,
                                         function(data) data$x2 > stats::median(evaluation$x2))
  oracle_target <- estimate_cspd(evaluation, "x1", grid, oracle, prediction)$cspd
  ale_on_grid <- stats::approx(ale$z, ale$ale, xout = grid, rule = 2)$y
  baseline_rows <- rbind(
    data.frame(method = "PDP", learner = NA_character_, curve_amplitude = s14_curve_amplitude(pdp$pdp),
               mean_isc = NA_real_, gamma = NA_real_, n_reference = n_evaluation,
               calibration_exceedance = NA_real_,
               shape_distance_to_oracle = s14_shape_distance(pdp$pdp, oracle_target)),
    data.frame(method = "conditional curve", learner = NA_character_, curve_amplitude = s14_curve_amplitude(conditional$conditional),
               mean_isc = NA_real_, gamma = NA_real_, n_reference = mean(conditional$n_context),
               calibration_exceedance = NA_real_,
               shape_distance_to_oracle = s14_shape_distance(conditional$conditional, oracle_target)),
    data.frame(method = "ALE", learner = NA_character_, curve_amplitude = s14_curve_amplitude(ale$ale),
               mean_isc = NA_real_, gamma = NA_real_, n_reference = sum(ale$n_bin),
               calibration_exceedance = NA_real_,
               shape_distance_to_oracle = s14_shape_distance(ale_on_grid, oracle_target)),
    data.frame(method = "subgroup PDP: low x2", learner = NA_character_, curve_amplitude = s14_curve_amplitude(low_subgroup$pdp),
               mean_isc = NA_real_, gamma = NA_real_, n_reference = unique(low_subgroup$n_subgroup),
               calibration_exceedance = NA_real_,
               shape_distance_to_oracle = s14_shape_distance(low_subgroup$pdp, oracle_target)),
    data.frame(method = "subgroup PDP: high x2", learner = NA_character_, curve_amplitude = s14_curve_amplitude(high_subgroup$pdp),
               mean_isc = NA_real_, gamma = NA_real_, n_reference = unique(high_subgroup$n_subgroup),
               calibration_exceedance = NA_real_,
               shape_distance_to_oracle = s14_shape_distance(high_subgroup$pdp, oracle_target))
  )
  cspd_rows <- do.call(rbind, lapply(names(learners), function(learner) {
    validity <- learners[[learner]]
    isc <- estimate_isc(evaluation, "x1", grid, validity)
    cspd <- estimate_cspd(evaluation, "x1", grid, validity, prediction)
    data.frame(method = "CSPD", learner = learner, curve_amplitude = s14_curve_amplitude(cspd$cspd),
               mean_isc = mean(isc$isc), gamma = unique(cspd$gamma), n_reference = unique(cspd$n_common),
               calibration_exceedance = calibration_exceedance[[learner]],
               shape_distance_to_oracle = s14_shape_distance(cspd$cspd, oracle_target))
  }))
  output <- rbind(baseline_rows, cspd_rows)
  transform(output, rho = rho, grid_half_width = grid_half_width,
            model_class = model_class, repetition = repetition,
            residual_density_identical = residual_density_identical)
}

run_s14 <- function(
  rhos = c(0.2, 0.5, 0.8),
  grid_half_widths = c(0.4, 0.75, 1.1),
  model_classes = c("linear", "tree", "gam", "neural_net", "random_forest", "boosted_tree"),
  repetitions = 100L,
  repetition_ids = seq_len(repetitions),
  n_training = 200L,
  n_calibration = 200L,
  n_evaluation = 400L,
  alpha = 0.05,
  k = 20L,
  workers = if (.Platform$OS.type == "windows") 1L else 8L,
  seed = 20260921L
) {
  design <- expand.grid(rho = rhos, grid_half_width = grid_half_widths,
                        model_class = model_classes, repetition = repetition_ids,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  evaluate_design <- function(index) {
    row <- design[index, ]
    s14_one(row$rho, row$grid_half_width, row$model_class, row$repetition,
            n_training, n_calibration, n_evaluation, alpha, k, seed + index)
  }
  indices <- seq_len(nrow(design))
  if (workers > 1L && .Platform$OS.type != "windows") {
    return(do.call(rbind, parallel::mclapply(indices, evaluate_design, mc.cores = workers)))
  }
  do.call(rbind, lapply(indices, evaluate_design))
}

summarize_s14 <- function(simulation) {
  simulation$learner[is.na(simulation$learner)] <- "not applicable"
  groups <- c("rho", "grid_half_width", "model_class", "method", "learner")
  measures <- c("curve_amplitude", "mean_isc", "gamma", "n_reference",
                "calibration_exceedance", "shape_distance_to_oracle")
  pieces <- split(simulation, interaction(simulation[groups], drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(pieces, function(part) {
    output <- part[1, groups, drop = FALSE]
    output$n_repetitions <- nrow(part)
    for (measure in measures) {
      values <- part[[measure]]
      values <- values[is.finite(values)]
      output[[measure]] <- if (length(values)) mean(values) else NA_real_
      output[[paste0(measure, "_mc_se")]] <- if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_
    }
    output
  }))
}

plot_s14 <- function(summary, path = file.path("results", "s14-experiment-matrix.png")) {
  cspd <- summary[summary$method == "CSPD", ]
  learners <- c("oracle", "split_conformal", "locally_scaled", "knn")
  cspd <- cspd[cspd$learner %in% learners, , drop = FALSE]
  learner_labels <- c(oracle = "Oracle", split_conformal = "Split residual",
                      locally_scaled = "Locally scaled", knn = "kNN")
  colours <- c("Oracle" = pdp_palette[["ink"]],
               "Split residual" = pdp_palette[["orange"]],
               "Locally scaled" = pdp_palette[["green"]],
               "kNN" = pdp_palette[["blue"]])
  shapes <- c("Oracle" = 16, "Split residual" = 17,
              "Locally scaled" = 8, "kNN" = 15)
  widest_grid <- max(cspd$grid_half_width)
  rho_values <- sort(unique(cspd$rho))
  offsets <- stats::setNames(seq(-0.025, 0.025, length.out = length(learners)), learners)
  pool_across_models <- function(data, measure) {
    parts <- split(data, interaction(data$rho, data$learner, drop = TRUE))
    do.call(rbind, lapply(parts, function(part) {
      means <- part[[measure]]
      sample_sizes <- part$n_repetitions
      standard_errors <- part[[paste0(measure, "_mc_se")]]
      standard_deviations <- standard_errors * sqrt(sample_sizes)
      pooled_mean <- stats::weighted.mean(means, sample_sizes)
      total_n <- sum(sample_sizes)
      pooled_variance <- if (total_n > 1) {
        sum((sample_sizes - 1) * standard_deviations^2 +
              sample_sizes * (means - pooled_mean)^2) / (total_n - 1)
      } else 0
      data.frame(rho = part$rho[[1]], learner = part$learner[[1]],
                 estimate = pooled_mean,
                 mc_se = sqrt(pooled_variance / total_n))
    }))
  }
  top_data <- do.call(rbind, lapply(c("gamma", "mean_isc"), function(measure) {
    pooled <- pool_across_models(cspd[cspd$grid_half_width == widest_grid, ], measure)
    pooled$measure <- measure
    pooled$x <- pooled$rho + offsets[pooled$learner]
    pooled$display_learner <- factor(unname(learner_labels[pooled$learner]),
                                     levels = unname(learner_labels[learners]))
    pooled$lower <- pmax(0, pooled$estimate - 2 * pooled$mc_se)
    pooled$upper <- pmin(1, pooled$estimate + 2 * pooled$mc_se)
    pooled
  }))
  top_panel <- function(measure, y_label, limits, breaks) {
    data <- top_data[top_data$measure == measure, ]
    ggplot2::ggplot(data, ggplot2::aes(x, estimate, color = display_learner,
                                      shape = display_learner,
                                      group = display_learner)) +
      ggplot2::geom_vline(xintercept = rho_values, color = "#E1E1E1",
                          linetype = "dotted", linewidth = 0.4) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 2.2, stroke = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                             width = 0.016, linewidth = 0.75) +
      ggplot2::scale_x_continuous(breaks = rho_values, labels = rho_values) +
      ggplot2::scale_y_continuous(breaks = breaks) +
      ggplot2::coord_cartesian(ylim = limits, clip = "off") +
      ggplot2::scale_color_manual(values = colours) +
      ggplot2::scale_shape_manual(values = shapes) +
      ggplot2::labs(x = expression(rho), y = y_label,
                    color = "Validity rule", shape = "Validity rule") +
      ggplot2::guides(color = "none", shape = "none") +
      theme_pdp() +
      ggplot2::theme(legend.position = "none")
  }
  p_gamma <- top_panel("gamma", expression(hat(gamma)(G)), c(0, 1.012), seq(0, 1, 0.2))
  p_isc <- top_panel("mean_isc", expression(widehat(ISC)^P(G)),
                     c(0.70, 1.006), seq(0.7, 1, 0.1))
  comparison <- cspd[cspd$learner %in% learners &
                       cspd$grid_half_width == widest_grid, ]
  comparison <- comparison[is.finite(comparison$gamma) &
                             is.finite(comparison$shape_distance_to_oracle), ]
  comparison$display_learner <- factor(
    unname(learner_labels[comparison$learner]),
    levels = unname(learner_labels[learners])
  )
  p_distance <- ggplot2::ggplot(
    comparison, ggplot2::aes(gamma, shape_distance_to_oracle,
                             color = display_learner, shape = display_learner)
  ) +
    ggplot2::geom_point(alpha = 0.38, size = 1.4) +
    ggplot2::scale_color_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(x = expression(hat(gamma)(G)),
                  y = "Centered distance to oracle CSPD",
                  color = "Validity rule", shape = "Validity rule") +
    ggplot2::guides(color = "none", shape = "none") +
    theme_pdp() +
    ggplot2::theme(legend.position = "none")
  model_data <- cspd[cspd$learner %in% learners &
                       cspd$grid_half_width == widest_grid &
                       cspd$rho == max(cspd$rho), ]
  model_labels <- c(
    boosted_tree = "Boosting", gam = "GAM", linear = "Linear",
    neural_net = "Neural net", random_forest = "Random forest", tree = "Tree"
  )
  model_data$display_model <- factor(
    unname(model_labels[model_data$model_class]),
    levels = unname(model_labels[c("linear", "tree", "gam", "neural_net",
                                  "random_forest", "boosted_tree")])
  )
  model_data$display_learner <- factor(
    unname(learner_labels[model_data$learner]),
    levels = unname(learner_labels[learners])
  )
  p_models <- ggplot2::ggplot(
    model_data,
    ggplot2::aes(display_model, shape_distance_to_oracle,
                 color = display_learner, shape = display_learner,
                 group = display_learner)
  ) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35),
                       linewidth = 0.45, alpha = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, shape_distance_to_oracle - 2 * shape_distance_to_oracle_mc_se),
        ymax = shape_distance_to_oracle + 2 * shape_distance_to_oracle_mc_se
      ),
      position = ggplot2::position_dodge(width = 0.35), width = 0.18,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.35),
                        size = 1.8) +
    ggplot2::scale_color_manual(values = colours) +
    ggplot2::scale_shape_manual(values = shapes) +
    ggplot2::labs(x = NULL, y = "Centered distance to oracle CSPD",
                  color = "Validity rule", shape = "Validity rule") +
    theme_pdp() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      legend.position = "bottom",
      legend.justification = "center"
    )
  save_pdp_panel_set(
    list(p_gamma, p_isc, p_distance, p_models), path,
    widths = 5.2, heights = 4.0
  )
}
