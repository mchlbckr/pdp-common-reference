# Lightweight screening of candidate real-data applications.
#
# This script is deliberately separate from the targets pipeline and the paper.
# It uses one fixed split, one random forest, K = 21, alpha = epsilon = 0.10,
# and at most 5,000 evaluation rows. Its purpose is triage, not final inference.

source(file.path("R", "support.R"))
source(file.path("R", "validity.R"))
source(file.path("R", "plot_theme.R"))

pilot_split <- function(n, seed) {
  set.seed(seed)
  index <- sample.int(n)
  n_training <- floor(0.60 * n)
  n_calibration <- floor(0.20 * n)
  list(
    training = index[seq_len(n_training)],
    calibration = index[n_training + seq_len(n_calibration)],
    evaluation = index[(n_training + n_calibration + 1L):n]
  )
}

pilot_cap_evaluation <- function(index, maximum = 5000L, seed = 20260916L) {
  if (length(index) <= maximum) return(index)
  set.seed(seed)
  sample(index, maximum)
}

pilot_standardizer <- function(training) {
  center <- vapply(training, mean, numeric(1))
  scale <- vapply(training, stats::sd, numeric(1))
  scale[!is.finite(scale) | scale <= 0] <- 1
  list(center = center, scale = scale)
}

pilot_standardize <- function(data, standardizer) {
  output <- sweep(as.matrix(data), 2, standardizer$center, "-")
  output <- sweep(output, 2, standardizer$scale, "/")
  as.data.frame(output, check.names = FALSE)
}

pilot_validity_rule <- function(training, calibration, focal, alpha = 0.10) {
  standardizer <- pilot_standardizer(training)
  fit <- suppressWarnings(fit_locally_scaled_validity(
    pilot_standardize(training, standardizer),
    pilot_standardize(calibration, standardizer),
    focal,
    alpha = alpha
  ))
  function(newdata) {
    suppressWarnings(predict_locally_scaled_validity(
      fit,
      pilot_standardize(newdata, standardizer)
    ))
  }
}

pilot_safe_distance <- function(first, second) {
  if (length(first) != length(second) || any(!is.finite(first)) ||
      any(!is.finite(second))) return(NA_real_)
  centered_curve_distance(first, second)
}

pilot_analyze_feature <- function(dataset, focal, training, calibration,
                                  evaluation, predict_function,
                                  observed_prediction_iqr,
                                  alpha = 0.10, epsilon = 0.10,
                                  grid_quantiles = c(0.30, 0.70), k = 21L) {
  is_valid <- pilot_validity_rule(training, calibration, focal, alpha)
  limits <- stats::quantile(
    training[[focal]], grid_quantiles, names = FALSE, type = 8
  )
  grid <- seq(limits[[1]], limits[[2]], length.out = k)
  validity <- validity_matrix(evaluation, focal, grid, is_valid)
  predictions <- prediction_matrix(evaluation, focal, grid, predict_function)
  pdp <- colMeans(predictions)
  common <- rowSums(validity) == ncol(validity)
  hard <- if (any(common)) colMeans(predictions[common, , drop = FALSE]) else
    rep(NA_real_, length(grid))
  pointwise <- vapply(seq_along(grid), function(j) {
    selected <- validity[, j]
    if (any(selected)) mean(predictions[selected, j]) else NA_real_
  }, numeric(1))
  epsilon_min <- minimum_relaxation_epsilon(validity)
  relaxed_fit <- if (epsilon + 2e-5 >= epsilon_min) {
    tryCatch(
      fit_relaxed_reference(validity, epsilon = epsilon, tolerance = 2e-5),
      error = function(error) NULL
    )
  } else {
    NULL
  }
  relaxed <- if (!is.null(relaxed_fit) && relaxed_fit$feasible) {
    as.numeric(crossprod(relaxed_fit$weights, predictions))
  } else {
    rep(NA_real_, length(grid))
  }
  prediction_iqr <- max(observed_prediction_iqr, sqrt(.Machine$double.eps))
  pd_amplitude <- diff(range(pdp - mean(pdp)))
  curves <- rbind(
    data.frame(z = grid, value = pdp, method = "Ordinary PD"),
    data.frame(z = grid, value = relaxed, method = "RCPD"),
    data.frame(z = grid, value = hard, method = "Hard CSPD"),
    data.frame(z = grid, value = pointwise, method = "PTPD")
  )
  curves$dataset <- dataset
  curves$focal <- focal
  curves$case <- paste(dataset, focal, sep = ": ")
  curves$centered <- ave(curves$value, curves$method, FUN = function(x) {
    if (all(is.finite(x))) x - mean(x) else x
  })
  curves$scaled_centered <- curves$centered / prediction_iqr

  support <- data.frame(
    z = grid,
    isc = colMeans(validity),
    dataset = dataset,
    focal = focal,
    case = paste(dataset, focal, sep = ": ")
  )
  summary <- data.frame(
    dataset = dataset,
    focal = focal,
    grid_lower = min(grid),
    grid_upper = max(grid),
    min_isc = min(colMeans(validity)),
    mean_isc = mean(validity),
    gamma = mean(common),
    epsilon_min = epsilon_min,
    rcpd_feasible = !is.null(relaxed_fit) && isTRUE(relaxed_fit$feasible),
    rcpd_ess = if (is.null(relaxed_fit)) NA_real_ else relaxed_fit$effective_sample_size,
    rcpd_kl = if (is.null(relaxed_fit)) NA_real_ else relaxed_fit$kl,
    pd_amplitude = pd_amplitude,
    prediction_iqr = observed_prediction_iqr,
    rcpd_distance = pilot_safe_distance(relaxed, pdp),
    hard_distance = pilot_safe_distance(hard, pdp),
    ptpd_distance = pilot_safe_distance(pointwise, pdp),
    rcpd_distance_over_pd_amplitude = pilot_safe_distance(relaxed, pdp) /
      max(pd_amplitude, sqrt(.Machine$double.eps)),
    rcpd_distance_over_prediction_iqr = pilot_safe_distance(relaxed, pdp) /
      prediction_iqr,
    pd_endpoint_contrast = tail(pdp, 1) - head(pdp, 1),
    rcpd_endpoint_contrast = if (all(is.finite(relaxed)))
      tail(relaxed, 1) - head(relaxed, 1) else NA_real_,
    hard_endpoint_contrast = if (all(is.finite(hard)))
      tail(hard, 1) - head(hard, 1) else NA_real_,
    ptpd_endpoint_contrast = if (all(is.finite(pointwise)))
      tail(pointwise, 1) - head(pointwise, 1) else NA_real_,
    stringsAsFactors = FALSE
  )
  list(summary = summary, curves = curves, support = support,
       validity = validity, predictions = predictions)
}

pilot_credit_data <- function() {
  path <- file.path(
    "data", "pilot", "credit", "default of credit card clients.csv"
  )
  raw <- utils::read.csv(path, skip = 1L, check.names = FALSE)
  names(raw)[names(raw) == "default payment next month"] <- "default"
  predictors <- raw[, setdiff(names(raw), c("ID", "default")), drop = FALSE]
  list(
    predictors = as.data.frame(lapply(predictors, as.numeric), check.names = FALSE),
    response = factor(ifelse(raw$default == 1, "yes", "no"),
                      levels = c("no", "yes"))
  )
}

pilot_california_data <- function() {
  path <- file.path(
    "data", "pilot", "california", "CaliforniaHousing", "cal_housing.data"
  )
  raw <- utils::read.csv(path, header = FALSE)
  names(raw) <- c(
    "Longitude", "Latitude", "HouseAge", "TotalRooms", "TotalBedrooms",
    "Population", "Households", "MedInc", "MedianHouseValue"
  )
  predictors <- transform(
    raw,
    AveRooms = TotalRooms / Households,
    AveBedrms = TotalBedrooms / Households,
    AveOccup = Population / Households
  )
  predictors <- predictors[, c(
    "MedInc", "HouseAge", "AveRooms", "AveBedrms", "Population",
    "AveOccup", "Latitude", "Longitude"
  )]
  list(
    predictors = predictors,
    response = raw$MedianHouseValue / 100000
  )
}

pilot_run_dataset <- function(dataset, data, focal_features, seed) {
  stopifnot(!anyNA(data$predictors), !anyNA(data$response))
  split <- pilot_split(nrow(data$predictors), seed)
  split$evaluation <- pilot_cap_evaluation(
    split$evaluation, maximum = 5000L, seed = seed + 1L
  )
  training <- data$predictors[split$training, , drop = FALSE]
  calibration <- data$predictors[split$calibration, , drop = FALSE]
  evaluation <- data$predictors[split$evaluation, , drop = FALSE]
  classification <- is.factor(data$response)
  set.seed(seed + 2L)
  model <- randomForest::randomForest(
    x = training,
    y = data$response[split$training],
    ntree = 100L,
    mtry = max(1L, floor(sqrt(ncol(training)))),
    nodesize = if (classification) 10L else 5L
  )
  predict_function <- if (classification) {
    function(newdata) as.numeric(stats::predict(
      model, newdata = newdata, type = "prob"
    )[, "yes"])
  } else {
    function(newdata) as.numeric(stats::predict(model, newdata = newdata))
  }
  observed_predictions <- predict_function(evaluation)
  response_evaluation <- data$response[split$evaluation]
  performance <- if (classification) {
    mean((as.numeric(response_evaluation == "yes") - observed_predictions)^2)
  } else {
    sqrt(mean((response_evaluation - observed_predictions)^2))
  }
  performance_name <- if (classification) "Brier" else "RMSE"
  prediction_iqr <- diff(stats::quantile(
    observed_predictions, c(0.25, 0.75), names = FALSE
  ))
  analyses <- lapply(focal_features, function(focal) {
    pilot_analyze_feature(
      dataset, focal, training, calibration, evaluation, predict_function,
      observed_prediction_iqr = prediction_iqr
    )
  })
  summaries <- do.call(rbind, lapply(analyses, `[[`, "summary"))
  summaries$n_training <- nrow(training)
  summaries$n_calibration <- nrow(calibration)
  summaries$n_evaluation <- nrow(evaluation)
  summaries$performance_name <- performance_name
  summaries$performance <- performance
  list(
    summary = summaries,
    curves = do.call(rbind, lapply(analyses, `[[`, "curves")),
    support = do.call(rbind, lapply(analyses, `[[`, "support")),
    model = model,
    analyses = analyses
  )
}

plot_pilot_curves <- function(curves, path) {
  curves$method <- factor(
    curves$method,
    levels = c("Ordinary PD", "RCPD", "Hard CSPD", "PTPD")
  )
  colors <- c(
    "Ordinary PD" = pdp_palette[["ink"]],
    "RCPD" = pdp_palette[["blue"]],
    "Hard CSPD" = pdp_palette[["orange"]],
    "PTPD" = pdp_palette[["green"]]
  )
  figure <- ggplot2::ggplot(
    curves,
    ggplot2::aes(z, scaled_centered, color = method, linetype = method)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = pdp_palette[["light_grey"]],
                        linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~ case, scales = "free_x", ncol = 2) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed", "dotdash", "dotted"), drop = FALSE
    ) +
    ggplot2::labs(
      x = "Focal-feature value",
      y = "Centered prediction / held-out prediction IQR",
      color = NULL,
      linetype = NULL
    ) +
    theme_pdp(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  save_pdp_plot(figure, path, width = 9.5, height = 7.0)
}

plot_pilot_support <- function(support, path) {
  figure <- ggplot2::ggplot(support, ggplot2::aes(z, isc)) +
    ggplot2::geom_hline(yintercept = 0.9, linetype = "dotted",
                        color = pdp_palette[["mid_grey"]]) +
    ggplot2::geom_line(color = pdp_palette[["green"]], linewidth = 0.8) +
    ggplot2::facet_wrap(~ case, scales = "free_x", ncol = 2) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(x = "Focal-feature value", y = "ISC under original reference") +
    theme_pdp(base_size = 10)
  save_pdp_plot(figure, path, width = 9.5, height = 7.0)
}

pilot_longitude_range_profile <- function(model, seed = 20260917L) {
  data <- pilot_california_data()
  split <- pilot_split(nrow(data$predictors), seed)
  split$evaluation <- pilot_cap_evaluation(
    split$evaluation, maximum = 5000L, seed = seed + 1L
  )
  training <- data$predictors[split$training, , drop = FALSE]
  calibration <- data$predictors[split$calibration, , drop = FALSE]
  evaluation <- data$predictors[split$evaluation, , drop = FALSE]
  predict_function <- function(newdata) {
    as.numeric(stats::predict(model, newdata = newdata))
  }
  observed_predictions <- predict_function(evaluation)
  prediction_iqr <- diff(stats::quantile(
    observed_predictions, c(0.25, 0.75), names = FALSE
  ))
  central_fractions <- c(0.10, 0.20, 0.30, 0.40)
  rows <- lapply(central_fractions, function(fraction) {
    result <- pilot_analyze_feature(
      "California housing", "Longitude", training, calibration, evaluation,
      predict_function, observed_prediction_iqr = prediction_iqr,
      grid_quantiles = c((1 - fraction) / 2, 1 - (1 - fraction) / 2)
    )
    transform(result$summary, central_fraction = fraction)
  })
  do.call(rbind, rows)
}

run_pilot_screen <- function() {
  credit <- pilot_run_dataset(
    "Credit default", pilot_credit_data(), c("LIMIT_BAL", "AGE"), 20260916L
  )
  california <- pilot_run_dataset(
    "California housing", pilot_california_data(),
    c("Longitude", "MedInc"), 20260917L
  )
  output <- list(
    summary = rbind(credit$summary, california$summary),
    curves = rbind(credit$curves, california$curves),
    support = rbind(credit$support, california$support),
    credit = credit,
    california = california
  )
  output$longitude_range_profile <- pilot_longitude_range_profile(
    california$model
  )
  dir.create(file.path("results", "pilot"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    output$summary,
    file.path("results", "pilot", "real-data-screen-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    output$longitude_range_profile,
    file.path("results", "pilot", "california-longitude-range-profile.csv"),
    row.names = FALSE
  )
  saveRDS(output, file.path("results", "pilot", "real-data-screen.rds"))
  plot_pilot_curves(
    output$curves,
    file.path("results", "pilot", "real-data-screen-curves.png")
  )
  plot_pilot_support(
    output$support,
    file.path("results", "pilot", "real-data-screen-support.png")
  )
  output
}

if (sys.nframe() == 0L) {
  screen <- run_pilot_screen()
  print(screen$summary, digits = 4)
}
