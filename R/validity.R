# Validity-domain estimators. The first estimator uses the standard
# split-conformal residual order statistic.

conformal_cutoff <- function(scores, alpha) {
  if (!is.numeric(scores) || length(scores) == 0L || any(!is.finite(scores))) {
    stop("scores must be a non-empty finite numeric vector.", call. = FALSE)
  }
  rank <- ceiling((length(scores) + 1) * (1 - alpha))
  if (rank > length(scores)) Inf else sort(scores, partial = rank)[rank]
}

fit_split_conformal_validity <- function(training_data, calibration_data, feature, alpha = 0.05) {
  assert_data_frame(training_data, "training_data")
  assert_data_frame(calibration_data, "calibration_data")
  assert_feature(training_data, feature)
  assert_feature(calibration_data, feature)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }

  predictors <- setdiff(names(training_data), feature)
  if (length(predictors) == 0L) stop("at least one conditioning feature is required.", call. = FALSE)
  formula <- stats::reformulate(predictors, response = feature)
  model <- stats::lm(formula, data = training_data)
  residuals <- abs(calibration_data[[feature]] - stats::predict(model, newdata = calibration_data))
  radius <- conformal_cutoff(residuals, alpha)

  structure(
    list(model = model, feature = feature, alpha = alpha, radius = radius,
         calibration_exceedance = mean(residuals > radius)),
    class = "split_conformal_validity"
  )
}

# Conditional Gaussian density-score validity.  The model and its residual scale
# are fit on training data; only the score threshold is calibrated separately.
fit_conditional_density_validity <- function(training_data, calibration_data, feature, alpha = 0.05) {
  assert_data_frame(training_data, "training_data")
  assert_data_frame(calibration_data, "calibration_data")
  assert_feature(training_data, feature)
  assert_feature(calibration_data, feature)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }

  predictors <- setdiff(names(training_data), feature)
  if (length(predictors) == 0L) stop("at least one conditioning feature is required.", call. = FALSE)
  formula <- stats::reformulate(predictors, response = feature)
  model <- stats::lm(formula, data = training_data)
  training_residuals <- stats::residuals(model)
  residual_scale <- sqrt(mean(training_residuals^2))
  if (!is.finite(residual_scale) || residual_scale <= 0) residual_scale <- .Machine$double.eps
  calibration_scores <- -stats::dnorm(
    calibration_data[[feature]],
    mean = stats::predict(model, newdata = calibration_data),
    sd = residual_scale,
    log = TRUE
  )
  threshold <- conformal_cutoff(calibration_scores, alpha)

  structure(
    list(model = model, feature = feature, alpha = alpha, residual_scale = residual_scale,
         threshold = threshold,
         calibration_exceedance = mean(calibration_scores > threshold)),
    class = "conditional_density_validity"
  )
}

predict_conditional_density_validity <- function(object, newdata) {
  if (!inherits(object, "conditional_density_validity")) stop("unsupported density validity model.", call. = FALSE)
  assert_data_frame(newdata, "newdata")
  if (!object$feature %in% names(newdata)) stop("newdata does not contain the focal feature.", call. = FALSE)
  scores <- -stats::dnorm(
    newdata[[object$feature]],
    mean = stats::predict(object$model, newdata = newdata),
    sd = object$residual_scale,
    log = TRUE
  )
  scores <= object$threshold
}

# Joint Gaussian compatibility score for numeric design matrices. Unlike a
# conditional focal-feature score, this rule also treats globally rare context
# vectors as unusual. A small ridge stabilizes one-hot encoded designs.
fit_joint_gaussian_validity <- function(training_data, calibration_data,
                                        alpha = 0.05, ridge = 1e-6) {
  assert_data_frame(training_data, "training_data")
  assert_data_frame(calibration_data, "calibration_data")
  if (!identical(names(training_data), names(calibration_data)) ||
      !all(vapply(training_data, is.numeric, logical(1))) ||
      !all(vapply(calibration_data, is.numeric, logical(1)))) {
    stop("training_data and calibration_data must have identical numeric columns.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge <= 0) {
    stop("ridge must be one positive finite number.", call. = FALSE)
  }
  center <- vapply(training_data, mean, numeric(1))
  covariance <- stats::cov(as.matrix(training_data))
  scale <- mean(diag(covariance)[diag(covariance) > 0])
  if (!is.finite(scale)) scale <- 1
  covariance <- covariance + diag(ridge * scale, ncol(training_data))
  precision <- chol2inv(chol(covariance))
  calibration_matrix <- sweep(as.matrix(calibration_data), 2, center, "-")
  scores <- rowSums((calibration_matrix %*% precision) * calibration_matrix)
  threshold <- conformal_cutoff(scores, alpha)
  structure(
    list(center = center, precision = precision, threshold = threshold,
         alpha = alpha, columns = names(training_data), ridge = ridge,
         calibration_exceedance = mean(scores > threshold)),
    class = "joint_gaussian_validity"
  )
}

predict_joint_gaussian_validity <- function(object, newdata) {
  if (!inherits(object, "joint_gaussian_validity")) {
    stop("unsupported joint Gaussian validity model.", call. = FALSE)
  }
  assert_data_frame(newdata, "newdata")
  if (!identical(names(newdata), object$columns)) {
    stop("newdata must have the fitted joint-score columns in the same order.", call. = FALSE)
  }
  centered <- sweep(as.matrix(newdata), 2, object$center, "-")
  scores <- rowSums((centered %*% object$precision) * centered)
  scores <= object$threshold
}

predict_validity <- function(object, newdata) {
  if (!inherits(object, "split_conformal_validity")) stop("unsupported validity model.", call. = FALSE)
  assert_data_frame(newdata, "newdata")
  if (!object$feature %in% names(newdata)) stop("newdata does not contain the focal feature.", call. = FALSE)
  center <- stats::predict(object$model, newdata = newdata)
  abs(newdata[[object$feature]] - center) <= object$radius
}

# Locally scaled split-conformal residual validity. This adapts interval widths
# to estimated heteroskedasticity while retaining marginal split-conformal
# calibration; it does not claim distribution-free conditional coverage.
fit_locally_scaled_validity <- function(training_data, calibration_data, feature, alpha = 0.05) {
  assert_data_frame(training_data, "training_data")
  assert_data_frame(calibration_data, "calibration_data")
  assert_feature(training_data, feature)
  assert_feature(calibration_data, feature)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between zero and one.", call. = FALSE)
  }
  predictors <- setdiff(names(training_data), feature)
  if (length(predictors) == 0L) stop("at least one conditioning feature is required.", call. = FALSE)
  formula <- stats::reformulate(predictors, response = feature)
  mean_model <- stats::lm(formula, data = training_data)
  scale_data <- training_data[predictors]
  scale_data$.log_abs_residual <- log(abs(stats::residuals(mean_model)) + sqrt(.Machine$double.eps))
  scale_formula <- stats::reformulate(predictors, response = ".log_abs_residual")
  scale_model <- stats::lm(scale_formula, data = scale_data)
  calibration_center <- stats::predict(mean_model, newdata = calibration_data)
  calibration_scale <- pmax(exp(stats::predict(scale_model, newdata = calibration_data)), sqrt(.Machine$double.eps))
  scores <- abs(calibration_data[[feature]] - calibration_center) / calibration_scale
  cutoff <- conformal_cutoff(scores, alpha)
  structure(
    list(mean_model = mean_model, scale_model = scale_model, feature = feature,
         alpha = alpha, cutoff = cutoff,
         calibration_exceedance = mean(scores > cutoff)),
    class = "locally_scaled_validity"
  )
}

predict_locally_scaled_validity <- function(object, newdata) {
  if (!inherits(object, "locally_scaled_validity")) stop("unsupported locally scaled validity model.", call. = FALSE)
  assert_data_frame(newdata, "newdata")
  if (!object$feature %in% names(newdata)) stop("newdata does not contain the focal feature.", call. = FALSE)
  center <- stats::predict(object$mean_model, newdata = newdata)
  scale <- pmax(exp(stats::predict(object$scale_model, newdata = newdata)), sqrt(.Machine$double.eps))
  abs(newdata[[object$feature]] - center) / scale <= object$cutoff
}

fit_knn_validity <- function(training_data, k = 25L, alpha = 0.05, calibration_data = NULL) {
  assert_data_frame(training_data, "training_data")
  if (k < 1L || k >= nrow(training_data)) stop("k must be between 1 and n-1.", call. = FALSE)
  center <- vapply(training_data, mean, numeric(1))
  scale <- vapply(training_data, stats::sd, numeric(1))
  scale[scale == 0] <- 1
  standardized <- sweep(sweep(as.matrix(training_data), 2, center, "-"), 2, scale, "/")
  if (is.null(calibration_data)) {
    distances <- as.matrix(stats::dist(standardized))
    diag(distances) <- Inf
    kth_distance <- apply(distances, 1, function(x) sort(x, partial = k)[k])
    threshold <- as.numeric(stats::quantile(kth_distance, probs = 1 - alpha, names = FALSE))
  } else {
    assert_data_frame(calibration_data, "calibration_data")
    if (!identical(names(calibration_data), names(training_data))) {
      stop("calibration_data must have the same columns as training_data.", call. = FALSE)
    }
    standardized_calibration <- sweep(sweep(as.matrix(calibration_data), 2, center, "-"), 2, scale, "/")
    kth_distance <- apply(standardized_calibration, 1, function(row) {
      values <- sqrt(rowSums((t(t(standardized) - row))^2))
      sort(values, partial = k)[k]
    })
    threshold <- conformal_cutoff(kth_distance, alpha)
  }
  structure(list(data = training_data, center = center, scale = scale, k = k,
                 threshold = threshold,
                 calibration_exceedance = mean(kth_distance > threshold)),
            class = "knn_validity")
}

predict_knn_validity <- function(object, newdata) {
  standardized_train <- sweep(sweep(as.matrix(object$data), 2, object$center, "-"), 2, object$scale, "/")
  standardized_new <- sweep(sweep(as.matrix(newdata), 2, object$center, "-"), 2, object$scale, "/")
  kth <- apply(standardized_new, 1, function(row) {
    distances <- sqrt(rowSums((t(t(standardized_train) - row))^2))
    sort(distances, partial = object$k)[object$k]
  })
  kth <= object$threshold
}
