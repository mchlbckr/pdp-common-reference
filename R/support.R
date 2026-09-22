# Support diagnostics and common-reference feature-effect estimators.
#
# This base-R implementation is the oracle implementation against which later
# optimized or package-facing code is tested.

assert_data_frame <- function(x, name) {
  if (!is.data.frame(x) || nrow(x) == 0L) {
    stop(sprintf("%s must be a non-empty data frame.", name), call. = FALSE)
  }
}

assert_feature <- function(data, feature) {
  if (!is.character(feature) || length(feature) != 1L || !feature %in% names(data)) {
    stop("feature must name one column of reference_data.", call. = FALSE)
  }
}

assert_grid <- function(grid) {
  if (!is.numeric(grid) || length(grid) == 0L || any(!is.finite(grid))) {
    stop("grid must be a non-empty vector of finite numeric values.", call. = FALSE)
  }
  unique(grid)
}

make_queries <- function(reference_data, feature, z) {
  query_data <- reference_data
  query_data[[feature]] <- z
  query_data
}

support_indicator <- function(reference_data, feature, z, is_valid) {
  indicator <- is_valid(make_queries(reference_data, feature, z))
  if (!is.logical(indicator) || length(indicator) != nrow(reference_data) || anyNA(indicator)) {
    stop("is_valid must return one non-missing logical value per query row.", call. = FALSE)
  }
  indicator
}

#' Estimate interventional support coverage (ISC) and its complement.
#'
#' @return A data frame with one row per grid value.
estimate_isc <- function(reference_data, feature, grid, is_valid) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)

  coverage <- vapply(
    grid,
    function(z) mean(support_indicator(reference_data, feature, z, is_valid)),
    numeric(1)
  )

  data.frame(z = grid, isc = coverage, per = 1 - coverage, stringsAsFactors = FALSE)
}

#' Estimate ordinary partial dependence on a supplied reference population.
estimate_pdp <- function(reference_data, feature, grid, predict_function) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)

  pdp <- vapply(grid, function(z) {
    predictions <- predict_function(make_queries(reference_data, feature, z))
    if (!is.numeric(predictions) || length(predictions) != nrow(reference_data) || any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per query row.", call. = FALSE)
    }
    mean(predictions)
  }, numeric(1))

  data.frame(z = grid, pdp = pdp, stringsAsFactors = FALSE)
}

#' Estimate partial dependence after pointwise validity trimming.
#'
#' Unlike CSPD and RCPD, the reference population is allowed to change with
#' `z`. The returned valid count and mass make that composition change
#' explicit rather than treating this curve as an ordinary PDP.
estimate_pointwise_trimmed_pdp <- function(reference_data, feature, grid,
                                           is_valid, predict_function) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)

  rows <- lapply(grid, function(z) {
    queries <- make_queries(reference_data, feature, z)
    selected <- is_valid(queries)
    if (!is.logical(selected) || length(selected) != nrow(reference_data) ||
        anyNA(selected)) {
      stop("is_valid must return one non-missing logical value per query row.",
           call. = FALSE)
    }
    n_valid <- sum(selected)
    if (n_valid == 0L) {
      return(c(pointwise_trimmed_pdp = NA_real_, n_valid = 0,
               valid_mass = 0))
    }
    predictions <- predict_function(queries[selected, , drop = FALSE])
    if (!is.numeric(predictions) || length(predictions) != n_valid ||
        any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per valid query row.",
           call. = FALSE)
    }
    c(pointwise_trimmed_pdp = mean(predictions), n_valid = n_valid,
      valid_mass = n_valid / nrow(reference_data))
  })
  output <- do.call(rbind, rows)
  data.frame(
    z = grid,
    pointwise_trimmed_pdp = output[, "pointwise_trimmed_pdp"],
    n_valid = as.integer(output[, "n_valid"]),
    valid_mass = output[, "valid_mass"],
    stringsAsFactors = FALSE
  )
}

assert_prediction_bounds <- function(prediction_bounds) {
  if (!is.numeric(prediction_bounds) || length(prediction_bounds) != 2L ||
      any(!is.finite(prediction_bounds)) ||
      prediction_bounds[[1]] >= prediction_bounds[[2]]) {
    stop("prediction_bounds must contain finite lower and upper values.",
         call. = FALSE)
  }
  unname(prediction_bounds)
}

#' Sharp empirical range-only bounds for ordinary partial dependence.
#'
#' Predictions are evaluated only on queries accepted by `is_valid`. Invalid
#' predictions are allowed to take any value in `prediction_bounds`.
estimate_range_pd_bounds <- function(reference_data, feature, grid, is_valid,
                                     predict_function, prediction_bounds,
                                     weights = NULL) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  prediction_bounds <- assert_prediction_bounds(prediction_bounds)
  n <- nrow(reference_data)
  if (is.null(weights)) weights <- rep(1 / n, n)
  if (!is.numeric(weights) || length(weights) != n ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("weights must be nonnegative, finite, sum to one, and match reference rows.",
         call. = FALSE)
  }

  rows <- lapply(grid, function(z) {
    queries <- make_queries(reference_data, feature, z)
    selected <- support_indicator(reference_data, feature, z, is_valid)
    predictions <- predict_function(queries[selected, , drop = FALSE])
    if (!is.numeric(predictions) || length(predictions) != sum(selected) ||
        any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per valid query row.",
           call. = FALSE)
    }
    if (length(predictions) &&
        any(predictions < prediction_bounds[[1]] - 1e-10 |
              predictions > prediction_bounds[[2]] + 1e-10)) {
      stop("valid-query predictions must lie within prediction_bounds.",
           call. = FALSE)
    }
    known_mean <- sum(weights[selected] * predictions)
    invalid_mass <- sum(weights[!selected])
    lower <- known_mean + prediction_bounds[[1]] * invalid_mass
    upper <- known_mean + prediction_bounds[[2]] * invalid_mass
    c(lower = lower, upper = upper, width = upper - lower,
      isc = sum(weights[selected]))
  })
  output <- do.call(rbind, rows)
  data.frame(z = grid, output, row.names = NULL, stringsAsFactors = FALSE)
}

euclidean_cross_distance <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (!is.numeric(x) || !is.numeric(y) || any(!is.finite(x)) ||
      any(!is.finite(y)) || ncol(x) != ncol(y)) {
    stop("the default context metric requires finite numeric context columns.",
         call. = FALSE)
  }
  if (ncol(x) == 0L) return(matrix(0, nrow(x), nrow(y)))
  squared <- outer(rowSums(x^2), rowSums(y^2), "+") -
    2 * tcrossprod(x, y)
  sqrt(pmax(squared, 0))
}

#' Sharp empirical Lipschitz-extension bounds for ordinary partial dependence.
#'
#' `distance_function` receives all contexts and the valid contexts and must
#' return their nonnegative cross-distance matrix. With `NULL`, raw Euclidean
#' distance over all numeric non-focal columns is used.
estimate_lipschitz_pd_bounds <- function(
  reference_data, feature, grid, is_valid, predict_function,
  lipschitz_constant, prediction_bounds, distance_function = NULL,
  tolerance = 1e-8, weights = NULL
) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  prediction_bounds <- assert_prediction_bounds(prediction_bounds)
  if (is.null(weights)) {
    weights <- rep(1 / nrow(reference_data), nrow(reference_data))
  }
  if (!is.numeric(weights) || length(weights) != nrow(reference_data) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("weights must be nonnegative, finite, sum to one, and match reference rows.",
         call. = FALSE)
  }
  if (!is.numeric(lipschitz_constant) || length(lipschitz_constant) != 1L ||
      !is.finite(lipschitz_constant) || lipschitz_constant < 0) {
    stop("lipschitz_constant must be one nonnegative finite number.",
         call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one nonnegative finite number.", call. = FALSE)
  }
  contexts <- reference_data[setdiff(names(reference_data), feature)]
  if (is.null(distance_function)) distance_function <- euclidean_cross_distance

  rows <- lapply(grid, function(z) {
    queries <- make_queries(reference_data, feature, z)
    selected <- support_indicator(reference_data, feature, z, is_valid)
    if (!any(selected)) {
      stop("each grid value needs at least one valid context for Lipschitz bounds.",
           call. = FALSE)
    }
    valid_predictions <- predict_function(queries[selected, , drop = FALSE])
    if (!is.numeric(valid_predictions) ||
        length(valid_predictions) != sum(selected) ||
        any(!is.finite(valid_predictions))) {
      stop("predict_function must return one finite numeric prediction per valid query row.",
           call. = FALSE)
    }
    if (any(valid_predictions < prediction_bounds[[1]] - tolerance |
            valid_predictions > prediction_bounds[[2]] + tolerance)) {
      stop("valid-query predictions must lie within prediction_bounds.",
           call. = FALSE)
    }
    distances <- distance_function(
      contexts, contexts[selected, , drop = FALSE]
    )
    if (!is.matrix(distances) || nrow(distances) != nrow(reference_data) ||
        ncol(distances) != sum(selected) || any(!is.finite(distances)) ||
        any(distances < 0)) {
      stop("distance_function must return a finite nonnegative n-by-n_valid matrix.",
           call. = FALSE)
    }
    valid_distances <- distances[selected, , drop = FALSE]
    pairwise_gaps <- abs(outer(valid_predictions, valid_predictions, "-"))
    if (any(pairwise_gaps > lipschitz_constant * valid_distances + tolerance)) {
      stop("valid-query predictions violate the declared Lipschitz constant.",
           call. = FALSE)
    }
    lower_candidates <- sweep(-lipschitz_constant * distances, 2,
                              valid_predictions, "+")
    upper_candidates <- sweep(lipschitz_constant * distances, 2,
                              valid_predictions, "+")
    lower_envelope <- pmax(prediction_bounds[[1]],
                           apply(lower_candidates, 1, max))
    upper_envelope <- pmin(prediction_bounds[[2]],
                           apply(upper_candidates, 1, min))
    lower <- sum(weights * lower_envelope)
    upper <- sum(weights * upper_envelope)
    c(lower = lower, upper = upper, width = upper - lower,
      isc = sum(weights[selected]))
  })
  output <- do.call(rbind, rows)
  data.frame(z = grid, output, row.names = NULL, stringsAsFactors = FALSE)
}

#' Estimate first-order accumulated local effects (ALE).
#'
#' This deliberately small implementation is used as a transparent baseline in
#' the simulations. `breaks` define the local intervals. Every observation must
#' lie in one interval; for observational data, quantile breaks are convenient.
#' The returned curve is evaluated at interval midpoints and centered with the
#' empirical interval frequencies, as in the usual first-order ALE definition.
estimate_ale <- function(reference_data, feature, breaks, predict_function) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  breaks <- assert_grid(breaks)
  if (length(breaks) < 2L || is.unsorted(breaks, strictly = TRUE)) {
    stop("breaks must contain at least two strictly increasing values.", call. = FALSE)
  }
  x <- reference_data[[feature]]
  if (any(x < breaks[[1]] | x > breaks[[length(breaks)]])) {
    stop("breaks must cover all observed feature values.", call. = FALSE)
  }

  bin <- findInterval(x, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  n_bins <- length(breaks) - 1L
  counts <- tabulate(bin, nbins = n_bins)
  if (any(counts == 0L)) {
    stop("each ALE interval must contain at least one observation.", call. = FALSE)
  }
  increments <- vapply(seq_len(n_bins), function(k) {
    rows <- bin == k
    upper <- make_queries(reference_data[rows, , drop = FALSE], feature, breaks[[k + 1L]])
    lower <- make_queries(reference_data[rows, , drop = FALSE], feature, breaks[[k]])
    mean(predict_function(upper) - predict_function(lower))
  }, numeric(1))
  uncentered <- cumsum(increments) - increments / 2
  centered <- uncentered - stats::weighted.mean(uncentered, counts)

  data.frame(
    z = (breaks[-1L] + breaks[-length(breaks)]) / 2,
    ale = centered,
    n_bin = counts,
    stringsAsFactors = FALSE
  )
}

#' Estimate a local conditional (M-plot-style) prediction curve.
#'
#' At every focal value, the reference population is the observed neighbourhood
#' of that value. This makes the context distribution vary with `z`; it is thus
#' deliberately not a fixed-reference interventional estimand.
estimate_conditional_curve <- function(reference_data, feature, grid, predict_function, bandwidth) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("bandwidth must be one positive finite number.", call. = FALSE)
  }
  x <- reference_data[[feature]]
  output <- lapply(grid, function(z) {
    rows <- abs(x - z) <= bandwidth
    if (!any(rows)) return(c(conditional = NA_real_, n_context = 0, context_mean = NA_real_))
    predictions <- predict_function(reference_data[rows, , drop = FALSE])
    c(conditional = mean(predictions), n_context = sum(rows), context_mean = NA_real_)
  })
  data.frame(
    z = grid,
    conditional = vapply(output, `[[`, numeric(1), "conditional"),
    n_context = vapply(output, `[[`, numeric(1), "n_context"),
    stringsAsFactors = FALSE
  )
}

#' Estimate ordinary PDP within a fixed, analyst-defined context subgroup.
#'
#' The subgroup is fixed over `grid`; this makes its population explicit, but
#' does not itself certify that all induced queries are valid.
estimate_subgroup_pdp <- function(reference_data, feature, grid, predict_function, subgroup_indicator) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  selected <- subgroup_indicator(reference_data)
  if (!is.logical(selected) || length(selected) != nrow(reference_data) || anyNA(selected) || !any(selected)) {
    stop("subgroup_indicator must select at least one reference row.", call. = FALSE)
  }
  estimate <- estimate_pdp(reference_data[selected, , drop = FALSE], feature, grid, predict_function)
  estimate$n_subgroup <- sum(selected)
  estimate$subgroup_mass <- mean(selected)
  estimate
}

#' Find contexts that support every intervention in a finite grid.
common_support_indicator <- function(reference_data, feature, grid, is_valid) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)

  Reduce(
    `&`,
    lapply(grid, function(z) support_indicator(reference_data, feature, z, is_valid))
  )
}

#' Evaluate a validity rule for every reference context and grid value.
#'
#' Rows index reference contexts and columns index interventions. This matrix is
#' the sufficient input for both hard and relaxed common-reference estimators.
validity_matrix <- function(reference_data, feature, grid, is_valid) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  matrix(
    unlist(lapply(grid, function(z) support_indicator(reference_data, feature, z, is_valid))),
    nrow = nrow(reference_data), ncol = length(grid),
    dimnames = list(NULL, format(grid, trim = TRUE))
  )
}

#' Critical intervention values for empirical interval validity.
#'
#' For row-wise closed intervals `[lower_i, upper_i]`, the validity pattern is
#' constant between consecutive endpoints. Its weighted coverage is upper
#' semicontinuous, so one midpoint from every open cell suffices for the exact
#' continuum minimum; endpoint constraints are redundant for relaxed validity.
interval_critical_grid <- function(lower, upper, interval) {
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) == 0L ||
      length(lower) != length(upper) || anyNA(lower) || anyNA(upper) ||
      any(lower > upper)) {
    stop("lower and upper must be equally long, non-empty numeric vectors with lower <= upper.",
         call. = FALSE)
  }
  if (!is.numeric(interval) || length(interval) != 2L ||
      any(!is.finite(interval)) || interval[[1]] >= interval[[2]]) {
    stop("interval must contain two finite, strictly increasing endpoints.",
         call. = FALSE)
  }
  inside <- c(lower[is.finite(lower)], upper[is.finite(upper)])
  inside <- inside[inside >= interval[[1]] & inside <= interval[[2]]]
  breakpoints <- sort(unique(c(interval, inside)))
  midpoints <- if (length(breakpoints) > 1L) {
    head(breakpoints, -1L) + diff(breakpoints) / 2
  } else {
    numeric()
  }
  sort(unique(midpoints))
}

#' Exact validity matrix for row-wise intervention intervals.
interval_validity_matrix <- function(lower, upper, grid) {
  grid <- assert_grid(grid)
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) == 0L ||
      length(lower) != length(upper) || anyNA(lower) || anyNA(upper) ||
      any(lower > upper)) {
    stop("lower and upper must be equally long, non-empty numeric vectors with lower <= upper.",
         call. = FALSE)
  }
  validity <- outer(lower, grid, `<=`) & outer(upper, grid, `>=`)
  dimnames(validity) <- list(NULL, format(grid, trim = TRUE))
  validity
}

#' Exact empirical semi-infinite KL projection for interval validity.
#'
#' The continuum of coverage constraints is reduced to the finitely many
#' validity patterns returned by `interval_critical_grid()`.
fit_interval_relaxed_reference <- function(
  lower, upper, interval, epsilon = 0.1, tolerance = 1e-7
) {
  critical_grid <- interval_critical_grid(lower, upper, interval)
  critical_validity <- interval_validity_matrix(lower, upper, critical_grid)
  fit <- fit_relaxed_reference(
    critical_validity, epsilon = epsilon, tolerance = tolerance
  )
  structure(list(
    weights = fit$weights,
    critical_grid = critical_grid,
    critical_validity = critical_validity,
    coverage = fit$coverage,
    min_coverage = min(fit$coverage),
    epsilon = epsilon,
    kl = fit$kl,
    effective_sample_size = fit$effective_sample_size,
    max_weight = fit$max_weight,
    fit = fit
  ), class = "interval_relaxed_reference")
}

#' Diagnostics common to every fixed-reference weighting rule.
reference_validity_diagnostics <- function(validity, weights = NULL,
                                           tolerance = 1e-8) {
  if (!is.matrix(validity) || !is.logical(validity) || nrow(validity) == 0L ||
      ncol(validity) == 0L || anyNA(validity)) {
    stop("validity must be a non-empty logical matrix without missing values.",
         call. = FALSE)
  }
  if (is.null(weights)) weights <- rep(1 / nrow(validity), nrow(validity))
  if (!is.numeric(weights) || length(weights) != nrow(validity) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > tolerance) {
    stop("weights must be nonnegative, finite, sum to one, and match validity rows.",
         call. = FALSE)
  }
  h <- unclass(validity) * 1
  coverage <- as.numeric(crossprod(weights, h))
  common <- rowSums(h) == ncol(h)
  list(
    coverage = coverage,
    min_coverage = min(coverage),
    mean_coverage = mean(coverage),
    simultaneous_valid_mass = sum(weights[common]),
    common_rows = sum(common)
  )
}

# KL projection subject to one expectation constraint E_Q score >= 1-epsilon.
# This supplies the average-validity and common-intersection comparators without
# altering the existing K-constraint RCPD implementation.
fit_single_constraint_reference <- function(validity, score, epsilon = 0.1,
                                            constraint = "single",
                                            tolerance = 1e-9) {
  if (!is.matrix(validity) || !is.logical(validity) || nrow(validity) == 0L ||
      ncol(validity) == 0L || anyNA(validity)) {
    stop("validity must be a non-empty logical matrix without missing values.",
         call. = FALSE)
  }
  if (!is.numeric(score) || length(score) != nrow(validity) ||
      any(!is.finite(score)) || any(score < 0) || any(score > 1)) {
    stop("score must be a finite vector in [0, 1] matching validity rows.",
         call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon < 0 || epsilon >= 1) {
    stop("epsilon must lie in [0, 1).", call. = FALSE)
  }
  n <- nrow(validity)
  target <- 1 - epsilon
  maximum <- max(score)
  epsilon_min <- 1 - maximum
  if (maximum + tolerance < target) {
    stop(sprintf(
      "%s constraint is infeasible; use epsilon >= %.6g.",
      constraint, epsilon_min
    ), call. = FALSE)
  }

  uniform <- rep(1 / n, n)
  achieved_uniform <- mean(score)
  if (achieved_uniform + tolerance >= target) {
    weights <- uniform
    multiplier <- 0
  } else if (abs(maximum - target) <= tolerance) {
    maximizers <- abs(score - maximum) <= tolerance
    weights <- maximizers / sum(maximizers)
    multiplier <- Inf
  } else {
    weighted_mean <- function(lambda) {
      eta <- lambda * score
      weights <- exp(eta - max(eta))
      weights <- weights / sum(weights)
      sum(weights * score)
    }
    upper <- 1
    while (weighted_mean(upper) < target && upper < 1e6) upper <- upper * 2
    multiplier <- stats::uniroot(
      function(lambda) weighted_mean(lambda) - target,
      interval = c(0, upper), tol = tolerance
    )$root
    eta <- multiplier * score
    weights <- exp(eta - max(eta))
    weights <- weights / sum(weights)
  }

  validity_diagnostics <- reference_validity_diagnostics(validity, weights)
  empirical_kl <- sum(ifelse(weights > 0, weights * log(weights * n), 0))
  pattern_diagnostics <- validity_pattern_diagnostics(validity, weights = weights)
  pattern_diagnostics$active_constraints <- as.integer(multiplier > tolerance)
  achieved <- sum(weights * score)
  structure(list(
    weights = weights,
    coverage = validity_diagnostics$coverage,
    min_coverage = validity_diagnostics$min_coverage,
    mean_coverage = validity_diagnostics$mean_coverage,
    simultaneous_valid_mass = validity_diagnostics$simultaneous_valid_mass,
    epsilon = epsilon,
    epsilon_min = max(0, epsilon_min),
    constraint = constraint,
    constraint_value = achieved,
    multiplier = multiplier,
    kl = empirical_kl,
    effective_sample_size = 1 / sum(weights^2),
    max_weight = max(weights),
    max_violation = max(0, target - achieved),
    feasible = achieved + tolerance >= target,
    pattern_diagnostics = pattern_diagnostics
  ), class = c("single_constraint_reference", "relaxed_reference"))
}

#' KL projection controlling mean pointwise validity over the grid.
fit_average_validity_reference <- function(validity, epsilon = 0.1,
                                           tolerance = 1e-9) {
  fit_single_constraint_reference(
    validity, rowMeans(validity), epsilon = epsilon,
    constraint = "average-validity", tolerance = tolerance
  )
}

#' KL projection controlling the mass valid simultaneously over the grid.
fit_intersection_reference <- function(validity, epsilon = 0.1,
                                       tolerance = 1e-9) {
  fit_single_constraint_reference(
    validity, as.numeric(rowSums(validity) == ncol(validity)),
    epsilon = epsilon, constraint = "intersection", tolerance = tolerance
  )
}

#' Evaluate a fixed set of reference weights over a PD grid.
estimate_weighted_pdp <- function(reference_data, feature, grid,
                                  predict_function, weights,
                                  value_name = "weighted_pdp") {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  if (!is.numeric(weights) || length(weights) != nrow(reference_data) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("weights must be nonnegative, finite, sum to one, and match reference rows.",
         call. = FALSE)
  }
  curve <- vapply(grid, function(z) {
    predictions <- predict_function(make_queries(reference_data, feature, z))
    if (!is.numeric(predictions) || length(predictions) != nrow(reference_data) ||
        any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per query row.",
           call. = FALSE)
    }
    sum(weights * predictions)
  }, numeric(1))
  output <- data.frame(z = grid, value = curve, stringsAsFactors = FALSE)
  names(output)[[2]] <- value_name
  output
}

#' Evaluate all induced-query predictions as a row-by-grid matrix.
prediction_matrix <- function(reference_data, feature, grid, predict_function) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  output <- vapply(grid, function(z) {
    predictions <- predict_function(make_queries(reference_data, feature, z))
    if (!is.numeric(predictions) || length(predictions) != nrow(reference_data) ||
        any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per query row.",
           call. = FALSE)
    }
    predictions
  }, numeric(nrow(reference_data)))
  if (is.null(dim(output))) output <- matrix(output, ncol = 1L)
  colnames(output) <- format(grid, trim = TRUE)
  output
}

#' Conditional evaluation-resampling stability for fixed learned objects.
#'
#' The validity and prediction matrices are treated as fixed row-level outputs.
#' Resampling therefore measures reference-sample and empirical-projection
#' variation, not training or validity-learning uncertainty.
bootstrap_reference_curve_stability <- function(
  validity, predictions, grid, epsilon = 0.1, repetitions = 300L,
  level = 0.95, seed = 20260929L,
  workers = if (.Platform$OS.type == "windows") 1L else 8L
) {
  if (!is.matrix(validity) || !is.logical(validity) || anyNA(validity)) {
    stop("validity must be a logical matrix without missing values.", call. = FALSE)
  }
  if (!is.matrix(predictions) || !is.numeric(predictions) ||
      any(!is.finite(predictions)) || !identical(dim(validity), dim(predictions))) {
    stop("predictions must be a finite numeric matrix matching validity.",
         call. = FALSE)
  }
  grid <- assert_grid(grid)
  if (length(grid) != ncol(validity)) {
    stop("grid length must match matrix columns.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("level must lie strictly between zero and one.", call. = FALSE)
  }
  n <- nrow(validity)
  curve_from_rows <- function(index) {
    h <- validity[index, , drop = FALSE]
    p <- predictions[index, , drop = FALSE]
    pd <- colMeans(p)
    common <- rowSums(h) == ncol(h)
    hard <- if (any(common)) colMeans(p[common, , drop = FALSE]) else
      rep(NA_real_, ncol(p))
    relaxed <- fit_relaxed_reference(h, epsilon = epsilon, tolerance = 2e-5)
    rcpd <- as.numeric(crossprod(relaxed$weights, p))
    list(
      curves = rbind(PDP = pd, `Hard CSPD` = hard, RCPD = rcpd),
      gamma = mean(common),
      ess = relaxed$effective_sample_size,
      min_coverage = relaxed$min_coverage,
      simultaneous_valid_mass = relaxed$simultaneous_valid_mass,
      normalized_hard_distance = if (all(is.finite(hard)))
        normalized_centered_curve_distance(hard, pd) else NA_real_,
      normalized_relaxed_distance = normalized_centered_curve_distance(rcpd, pd)
    )
  }
  original <- curve_from_rows(seq_len(n))
  set.seed(seed)
  indices <- replicate(repetitions, sample.int(n, n, replace = TRUE),
                       simplify = FALSE)
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(indices, curve_from_rows, mc.cores = workers)
  } else {
    lapply(indices, curve_from_rows)
  }
  method_names <- rownames(original$curves)
  center_rows <- function(matrix) matrix - rowMeans(matrix)
  original_centered <- center_rows(original$curves)
  curve_array <- array(
    unlist(lapply(pieces, `[[`, "curves")),
    dim = c(length(method_names), length(grid), repetitions),
    dimnames = list(method_names, NULL, NULL)
  )
  bands <- do.call(rbind, lapply(seq_along(method_names), function(index) {
    draws <- t(vapply(seq_len(repetitions), function(b) {
      values <- curve_array[index, , b]
      values - mean(values)
    }, numeric(length(grid))))
    point <- original_centered[index, ]
    maximum_deviation <- apply(abs(sweep(draws, 2, point, "-")), 1, max,
                               na.rm = TRUE)
    critical <- stats::quantile(maximum_deviation, level, na.rm = TRUE,
                                names = FALSE)
    data.frame(
      z = grid, method = method_names[[index]], estimate = point,
      median = apply(draws, 2, stats::median, na.rm = TRUE),
      lower = point - critical, upper = point + critical,
      level = level, stringsAsFactors = FALSE
    )
  }))
  differences <- do.call(rbind, lapply(c("Hard CSPD", "RCPD"), function(method) {
    index <- match(method, method_names)
    pd_index <- match("PDP", method_names)
    point <- original_centered[index, ] - original_centered[pd_index, ]
    draws <- t(vapply(seq_len(repetitions), function(b) {
      first <- curve_array[index, , b]
      second <- curve_array[pd_index, , b]
      (first - mean(first)) - (second - mean(second))
    }, numeric(length(grid))))
    maximum_deviation <- apply(abs(sweep(draws, 2, point, "-")), 1, max,
                               na.rm = TRUE)
    critical <- stats::quantile(maximum_deviation, level, na.rm = TRUE,
                                names = FALSE)
    data.frame(
      z = grid, contrast = paste(method, "minus PDP"), estimate = point,
      lower = point - critical, upper = point + critical,
      level = level, stringsAsFactors = FALSE
    )
  }))
  metrics <- data.frame(
    repetition = seq_len(repetitions),
    gamma = vapply(pieces, `[[`, numeric(1), "gamma"),
    ess = vapply(pieces, `[[`, numeric(1), "ess"),
    min_coverage = vapply(pieces, `[[`, numeric(1), "min_coverage"),
    simultaneous_valid_mass = vapply(
      pieces, `[[`, numeric(1), "simultaneous_valid_mass"
    ),
    normalized_hard_distance = vapply(
      pieces, `[[`, numeric(1), "normalized_hard_distance"
    ),
    normalized_relaxed_distance = vapply(
      pieces, `[[`, numeric(1), "normalized_relaxed_distance"
    )
  )
  list(bands = bands, differences = differences, metrics = metrics,
       original = original, repetitions = repetitions, level = level)
}

#' Solve the grid-wise KL projection by adaptive constraint generation.
#'
#' The projection starts from a small grid, evaluates its achieved coverage on
#' a dense validation grid, and repeatedly adds the worst violated point. This
#' approximates a continuum-indexed minimum-coverage constraint without making
#' the initial display grid part of the estimand.
fit_adaptive_relaxed_reference <- function(
  reference_data, feature, interval, is_valid, epsilon = 0.1,
  initial_grid = NULL, validation_grid = NULL,
  tolerance = 2e-5, max_iterations = 100L
) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  if (!is.numeric(interval) || length(interval) != 2L ||
      any(!is.finite(interval)) || interval[[1]] >= interval[[2]]) {
    stop("interval must contain two finite, strictly increasing endpoints.",
         call. = FALSE)
  }
  if (is.null(validation_grid)) {
    validation_grid <- seq(interval[[1]], interval[[2]], length.out = 201L)
  }
  validation_grid <- sort(assert_grid(validation_grid))
  if (min(validation_grid) < interval[[1]] ||
      max(validation_grid) > interval[[2]]) {
    stop("validation_grid must lie inside interval.", call. = FALSE)
  }
  if (is.null(initial_grid)) {
    initial_grid <- c(interval[[1]], mean(interval), interval[[2]])
  }
  active_grid <- sort(unique(c(interval, assert_grid(initial_grid))))
  if (min(active_grid) < interval[[1]] || max(active_grid) > interval[[2]]) {
    stop("initial_grid must lie inside interval.", call. = FALSE)
  }
  full_validity <- validity_matrix(
    reference_data, feature, validation_grid, is_valid
  )
  converged <- FALSE
  fit <- NULL
  validation_coverage <- NULL
  for (iteration in seq_len(max_iterations)) {
    active_validity <- validity_matrix(
      reference_data, feature, active_grid, is_valid
    )
    fit <- fit_relaxed_reference(
      active_validity, epsilon = epsilon, tolerance = tolerance
    )
    validation_coverage <- as.numeric(crossprod(
      fit$weights, unclass(full_validity) * 1
    ))
    worst <- which.min(validation_coverage)
    if (validation_coverage[[worst]] >= 1 - epsilon - tolerance) {
      converged <- TRUE
      break
    }
    candidate <- validation_grid[[worst]]
    if (candidate %in% active_grid) break
    active_grid <- sort(c(active_grid, candidate))
  }
  diagnostics <- reference_validity_diagnostics(full_validity, fit$weights)
  structure(list(
    weights = fit$weights,
    active_grid = active_grid,
    validation_grid = validation_grid,
    validation_coverage = validation_coverage,
    min_validation_coverage = min(validation_coverage),
    mean_validation_coverage = mean(validation_coverage),
    simultaneous_valid_mass = diagnostics$simultaneous_valid_mass,
    iterations = iteration,
    converged = converged,
    epsilon = epsilon,
    kl = fit$kl,
    effective_sample_size = fit$effective_sample_size,
    max_weight = max(fit$weights),
    fit = fit
  ), class = "adaptive_relaxed_reference")
}

#' Solve empirical continuum constraints by exact separation.
#'
#' This exchange algorithm is exact for row-wise interval validity. Its
#' separation step searches every distinct validity pattern over the declared
#' interval, so finite termination follows once no critical pattern violates
#' the requested coverage.
fit_adaptive_interval_relaxed_reference <- function(
  lower, upper, interval, epsilon = 0.1, initial_grid = NULL,
  tolerance = 2e-5, max_iterations = NULL
) {
  critical_grid <- interval_critical_grid(lower, upper, interval)
  full_validity <- interval_validity_matrix(lower, upper, critical_grid)
  if (is.null(initial_grid)) {
    initial_grid <- c(interval[[1]], mean(interval), interval[[2]])
  }
  active_grid <- sort(unique(assert_grid(initial_grid)))
  if (min(active_grid) < interval[[1]] || max(active_grid) > interval[[2]]) {
    stop("initial_grid must lie inside interval.", call. = FALSE)
  }
  if (is.null(max_iterations)) max_iterations <- length(critical_grid)
  if (!is.numeric(max_iterations) || length(max_iterations) != 1L ||
      !is.finite(max_iterations) || max_iterations < 1L) {
    stop("max_iterations must be one positive finite number.", call. = FALSE)
  }

  converged <- FALSE
  fit <- NULL
  critical_coverage <- NULL
  for (iteration in seq_len(as.integer(max_iterations))) {
    active_validity <- interval_validity_matrix(lower, upper, active_grid)
    fit <- fit_relaxed_reference(
      active_validity, epsilon = epsilon, tolerance = tolerance
    )
    critical_coverage <- as.numeric(crossprod(
      fit$weights, unclass(full_validity) * 1
    ))
    worst <- which.min(critical_coverage)
    if (critical_coverage[[worst]] >= 1 - epsilon - tolerance) {
      converged <- TRUE
      break
    }
    candidate <- critical_grid[[worst]]
    if (candidate %in% active_grid) break
    active_grid <- sort(c(active_grid, candidate))
  }
  diagnostics <- reference_validity_diagnostics(full_validity, fit$weights)
  structure(list(
    weights = fit$weights,
    active_grid = active_grid,
    critical_grid = critical_grid,
    critical_coverage = critical_coverage,
    min_critical_coverage = min(critical_coverage),
    mean_critical_coverage = mean(critical_coverage),
    simultaneous_valid_mass = diagnostics$simultaneous_valid_mass,
    iterations = iteration,
    converged = converged,
    epsilon = epsilon,
    kl = fit$kl,
    effective_sample_size = fit$effective_sample_size,
    max_weight = max(fit$weights),
    fit = fit
  ), class = "adaptive_interval_relaxed_reference")
}

#' Smallest feasible relaxed-validity tolerance on an empirical pattern simplex.
#'
#' Solves min_q max_k sum_i q_i (1 - H_ik) over probability weights q. The
#' calculation is performed after collapsing duplicate validity patterns.
minimum_relaxation_epsilon <- function(validity) {
  if (!is.matrix(validity) || !is.logical(validity) || nrow(validity) == 0L ||
      ncol(validity) == 0L || anyNA(validity)) {
    stop("validity must be a non-empty logical matrix without missing values.", call. = FALSE)
  }
  patterns <- unique(unclass(validity) * 1)
  n_patterns <- nrow(patterns)
  n_grid <- ncol(patterns)
  objective <- c(rep(0, n_patterns), 1)
  constraints <- cbind(1 - t(patterns), -1)
  fit <- lpSolve::lp(
    direction = "min", objective.in = objective,
    const.mat = rbind(constraints, c(rep(1, n_patterns), 0)),
    const.dir = c(rep("<=", n_grid), "="),
    const.rhs = c(rep(0, n_grid), 1),
    all.int = FALSE
  )
  if (fit$status != 0L) {
    stop("could not solve the empirical minimum-epsilon linear program.", call. = FALSE)
  }
  max(0, min(1, fit$objval))
}

#' Summarize empirical validity patterns and a relaxed reference fit.
validity_pattern_diagnostics <- function(validity, weights = NULL,
                                         multipliers = NULL, tolerance = 1e-7) {
  if (!is.matrix(validity) || !is.logical(validity) || nrow(validity) == 0L ||
      ncol(validity) == 0L || anyNA(validity)) {
    stop("validity must be a non-empty logical matrix without missing values.", call. = FALSE)
  }
  if (is.null(weights)) weights <- rep(1 / nrow(validity), nrow(validity))
  if (!is.numeric(weights) || length(weights) != nrow(validity) ||
      any(!is.finite(weights)) || any(weights < 0) || abs(sum(weights) - 1) > tolerance) {
    stop("weights must be nonnegative, finite, sum to one, and match validity rows.", call. = FALSE)
  }
  keys <- apply(validity, 1, function(row) paste0(as.integer(row), collapse = ""))
  cell_sizes <- as.numeric(table(keys))
  active_constraints <- if (is.null(multipliers)) {
    NA_integer_
  } else {
    if (!is.numeric(multipliers) || length(multipliers) != ncol(validity)) {
      stop("multipliers must be numeric and match validity columns.", call. = FALSE)
    }
    sum(is.finite(multipliers) & multipliers > tolerance)
  }
  active_affine_dimension <- if (is.null(multipliers)) {
    NA_integer_
  } else {
    active <- is.finite(multipliers) & multipliers > tolerance
    if (!any(active)) {
      0L
    } else {
      active_patterns <- unique((unclass(validity[, active, drop = FALSE]) * 1))
      if (nrow(active_patterns) <= 1L) 0L else
        qr(sweep(active_patterns[-1, , drop = FALSE], 2L,
                 active_patterns[1, ], "-"), tol = tolerance)$rank
    }
  }
  active_licq <- if (is.na(active_constraints)) {
    NA
  } else {
    active_affine_dimension == active_constraints
  }
  data.frame(
    observed_patterns = length(cell_sizes),
    median_pattern_size = stats::median(cell_sizes),
    min_pattern_size = min(cell_sizes),
    singleton_pattern_fraction = mean(cell_sizes == 1L),
    active_constraints = active_constraints,
    active_affine_dimension = active_affine_dimension,
    active_licq = active_licq,
    max_weight = max(weights),
    epsilon_min = minimum_relaxation_epsilon(validity),
    stringsAsFactors = FALSE
  )
}

log_mean_exp <- function(x) {
  anchor <- max(x)
  anchor + log(mean(exp(x - anchor)))
}

#' Certify a finite family of RCPD tilts on an independent audit split.
#'
#' Candidate multipliers are learned only from `projection_validity`. Conditional
#' on that split, simultaneous one-sided bounds over candidates and grid points
#' permit audit-data-dependent selection without invalidating the certificate.
certify_calibrated_rcpd <- function(
  projection_validity, audit_validity,
  fitted_epsilons = c(0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10),
  target_epsilon = 0.10, delta = 0.10,
  bound = c("empirical_bernstein", "hoeffding"),
  tolerance = 2e-5
) {
  bound <- match.arg(bound)
  valid_matrix <- function(x) {
    is.matrix(x) && is.logical(x) && nrow(x) > 0L && ncol(x) > 0L &&
      !anyNA(x)
  }
  if (!valid_matrix(projection_validity) || !valid_matrix(audit_validity)) {
    stop("projection_validity and audit_validity must be non-empty logical matrices.",
         call. = FALSE)
  }
  if (ncol(projection_validity) != ncol(audit_validity)) {
    stop("projection and audit validity matrices must have the same columns.",
         call. = FALSE)
  }
  if (!is.numeric(target_epsilon) || length(target_epsilon) != 1L ||
      !is.finite(target_epsilon) || target_epsilon <= 0 ||
      target_epsilon >= 1) {
    stop("target_epsilon must lie strictly between zero and one.",
         call. = FALSE)
  }
  fitted_epsilons <- sort(unique(fitted_epsilons))
  if (!is.numeric(fitted_epsilons) || length(fitted_epsilons) == 0L ||
      any(!is.finite(fitted_epsilons)) || any(fitted_epsilons <= 0) ||
      any(fitted_epsilons > target_epsilon)) {
    stop("fitted_epsilons must be positive and no larger than target_epsilon.",
         call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta <= 0 || delta >= 1) {
    stop("delta must lie strictly between zero and one.", call. = FALSE)
  }
  if (bound == "empirical_bernstein" && nrow(audit_validity) < 2L) {
    stop("empirical Bernstein certification needs at least two audit rows.",
         call. = FALSE)
  }

  candidate_count <- length(fitted_epsilons)
  grid_size <- ncol(projection_validity)
  audit_size <- nrow(audit_validity)
  target <- 1 - target_epsilon
  audit_numeric <- unclass(audit_validity) * 1
  family_size <- candidate_count * grid_size
  fits <- vector("list", candidate_count)
  candidate_rows <- vector("list", candidate_count)
  constraint_rows <- vector("list", candidate_count)

  for (candidate in seq_along(fitted_epsilons)) {
    fitted_epsilon <- fitted_epsilons[[candidate]]
    fit <- tryCatch(
      fit_relaxed_reference(
        projection_validity, epsilon = fitted_epsilon,
        tolerance = tolerance
      ),
      error = identity
    )
    if (inherits(fit, "error") || !fit$feasible ||
        any(!is.finite(fit$multipliers))) {
      candidate_rows[[candidate]] <- data.frame(
        candidate = candidate, fitted_epsilon = fitted_epsilon,
        status = "projection-infeasible", certified = FALSE,
        projection_kl = NA_real_, projection_ess = NA_real_,
        audit_min_coverage = NA_real_, min_lcb_margin = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    multipliers <- pmax(fit$multipliers, 0)
    eta <- as.numeric(audit_numeric %*% multipliers)
    tilt <- exp(eta - sum(multipliers))
    moment_values <- tilt * sweep(audit_numeric, 2, target, "-")
    moment_means <- colMeans(moment_values)
    range_widths <- (1 - target) + target * exp(-multipliers)
    empirical_variance <- apply(moment_values, 2, stats::var)
    if (bound == "hoeffding") {
      log_term <- log(family_size / delta)
      radii <- range_widths * sqrt(log_term / (2 * audit_size))
    } else {
      log_term <- log(2 * family_size / delta)
      radii <- sqrt(2 * empirical_variance * log_term / audit_size) +
        7 * range_widths * log_term / (3 * (audit_size - 1))
    }
    lcb_margin <- moment_means - radii
    normalized_weights <- tilt / sum(tilt)
    audit_coverage <- as.numeric(crossprod(normalized_weights, audit_numeric))
    certified <- all(lcb_margin >= 0)

    fits[[candidate]] <- fit
    constraint_rows[[candidate]] <- data.frame(
      candidate = candidate,
      fitted_epsilon = fitted_epsilon,
      grid_index = seq_len(grid_size),
      multiplier = multipliers,
      audit_coverage = audit_coverage,
      audit_moment = moment_means,
      empirical_variance = empirical_variance,
      range_width = range_widths,
      radius = radii,
      lcb_margin = lcb_margin,
      stringsAsFactors = FALSE
    )
    candidate_rows[[candidate]] <- data.frame(
      candidate = candidate, fitted_epsilon = fitted_epsilon,
      status = "audited", certified = certified,
      projection_kl = fit$kl,
      projection_ess = fit$effective_sample_size,
      audit_min_coverage = min(audit_coverage),
      min_lcb_margin = min(lcb_margin),
      stringsAsFactors = FALSE
    )
  }

  candidates <- do.call(rbind, candidate_rows)
  constraints <- do.call(
    rbind, constraint_rows[!vapply(constraint_rows, is.null, logical(1))]
  )
  eligible <- which(candidates$certified)
  selected <- NULL
  if (length(eligible)) {
    ordering <- order(
      candidates$projection_kl[eligible],
      -candidates$fitted_epsilon[eligible]
    )
    selected_index <- candidates$candidate[eligible[ordering[[1]]]]
    selected <- list(
      candidate = selected_index,
      fitted_epsilon = fitted_epsilons[[selected_index]],
      fit = fits[[selected_index]],
      multipliers = pmax(fits[[selected_index]]$multipliers, 0),
      audit = candidates[candidates$candidate == selected_index, , drop = FALSE]
    )
  }

  structure(list(
    candidates = candidates,
    constraints = constraints,
    selected = selected,
    bound = bound,
    target_epsilon = target_epsilon,
    target_coverage = target,
    delta = delta,
    candidate_count = candidate_count,
    grid_size = grid_size,
    audit_size = audit_size
  ), class = "calibrated_rcpd_certificate")
}

evaluate_rcpd_tilt <- function(validity, multipliers) {
  weights <- apply_relaxed_reference_tilt(validity, multipliers)
  coverage <- as.numeric(crossprod(weights, unclass(validity) * 1))
  list(
    coverage = coverage,
    min_coverage = min(coverage),
    mean_coverage = mean(coverage),
    effective_sample_fraction = 1 / sum(weights^2) / nrow(validity)
  )
}

#' Apply fitted relaxed-reference multipliers to independent contexts.
#'
#' This is useful for an out-of-sample audit: the multipliers are learned on a
#' projection sample, while the returned normalized exponential-tilt weights
#' are computed from a new validity matrix.
apply_relaxed_reference_tilt <- function(validity, multipliers) {
  if (!is.matrix(validity) || !is.logical(validity) || nrow(validity) == 0L ||
      ncol(validity) == 0L || anyNA(validity)) {
    stop("validity must be a non-empty logical matrix without missing values.",
         call. = FALSE)
  }
  if (!is.numeric(multipliers) || length(multipliers) != ncol(validity) ||
      any(!is.finite(multipliers)) || any(multipliers < -1e-10)) {
    stop("multipliers must be finite, nonnegative, and match validity columns.",
         call. = FALSE)
  }
  # L-BFGS-B may return a boundary coordinate a few ulps below zero.
  multipliers <- pmax(multipliers, 0)
  eta <- as.numeric((unclass(validity) * 1) %*% multipliers)
  weights <- exp(eta - max(eta))
  weights / sum(weights)
}

#' Empirical KL projection under relaxed grid-wise validity constraints.
#'
#' Finds weights q_i closest to uniform empirical weights in KL divergence,
#' subject to sum_i q_i H_ik >= 1 - epsilon for every grid point k. The dual
#' solution has q_i proportional to exp(sum_k lambda_k H_ik).
fit_relaxed_reference <- function(validity, epsilon = 0.1, tolerance = 1e-7,
                                  max_multiplier = 50, max_iterations = 5000L) {
  if (!is.matrix(validity) || !is.logical(validity) || nrow(validity) == 0L ||
      ncol(validity) == 0L || anyNA(validity)) {
    stop("validity must be a non-empty logical matrix without missing values.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon < 0 || epsilon >= 1) {
    stop("epsilon must lie in [0, 1).", call. = FALSE)
  }
  n <- nrow(validity)
  k <- ncol(validity)
  h <- unclass(validity) * 1
  epsilon_min <- minimum_relaxation_epsilon(validity)
  if (epsilon + tolerance < epsilon_min) {
    stop(sprintf(
      "epsilon = %.6g is infeasible for the observed validity patterns; use epsilon >= %.6g.",
      epsilon, epsilon_min
    ), call. = FALSE)
  }

  if (epsilon == 0) {
    common <- rowSums(h) == k
    if (!any(common)) {
      stop("the hard reference is infeasible because no context is valid over the full grid.", call. = FALSE)
    }
    weights <- common / sum(common)
    coverage <- as.numeric(crossprod(weights, h))
    diagnostics <- validity_pattern_diagnostics(validity, weights = weights)
    validity_diagnostics <- reference_validity_diagnostics(validity, weights)
    return(structure(list(
      weights = weights, coverage = coverage, epsilon = epsilon,
      multipliers = rep(Inf, k), kl = -log(mean(common)),
      effective_sample_size = sum(common), max_violation = 0,
      convergence = 0L, feasible = TRUE, epsilon_min = epsilon_min,
      min_coverage = validity_diagnostics$min_coverage,
      mean_coverage = validity_diagnostics$mean_coverage,
      simultaneous_valid_mass = validity_diagnostics$simultaneous_valid_mass,
      binding_constraints = k, strict_complementarity = NA,
      pattern_diagnostics = diagnostics
    ), class = "relaxed_reference"))
  }

  target <- 1 - epsilon
  objective <- function(lambda) {
    eta <- as.numeric(h %*% lambda)
    log_mean_exp(eta) - target * sum(lambda)
  }
  gradient <- function(lambda) {
    eta <- as.numeric(h %*% lambda)
    weights <- exp(eta - max(eta))
    weights <- weights / sum(weights)
    as.numeric(crossprod(weights, h)) - target
  }
  fit <- stats::optim(
    par = rep(0, k), fn = objective, gr = gradient, method = "L-BFGS-B",
    lower = rep(0, k), upper = rep(max_multiplier, k),
    control = list(maxit = max_iterations, factr = 1e7)
  )
  eta <- as.numeric(h %*% fit$par)
  weights <- exp(eta - max(eta))
  weights <- weights / sum(weights)
  coverage <- as.numeric(crossprod(weights, h))
  max_violation <- max(0, target - min(coverage))
  empirical_kl <- sum(weights * log(weights * n))
  effective_sample_size <- 1 / sum(weights^2)
  binding <- abs(coverage - target) <= sqrt(tolerance)
  multiplier_active <- fit$par > sqrt(tolerance)
  strict_complementarity <- all(binding == multiplier_active)
  diagnostics <- validity_pattern_diagnostics(
    validity, weights = weights, multipliers = fit$par, tolerance = tolerance
  )
  validity_diagnostics <- reference_validity_diagnostics(validity, weights)

  structure(list(
    weights = weights, coverage = coverage, epsilon = epsilon,
    multipliers = fit$par, kl = empirical_kl,
    effective_sample_size = effective_sample_size,
    max_violation = max_violation, convergence = fit$convergence,
    feasible = fit$convergence == 0L && max_violation <= tolerance,
    epsilon_min = epsilon_min, binding_constraints = sum(binding),
    min_coverage = validity_diagnostics$min_coverage,
    mean_coverage = validity_diagnostics$mean_coverage,
    simultaneous_valid_mass = validity_diagnostics$simultaneous_valid_mass,
    strict_complementarity = strict_complementarity,
    pattern_diagnostics = diagnostics
  ), class = "relaxed_reference")
}

#' Estimate relaxed common-reference partial dependence (RC-PD).
estimate_relaxed_cspd <- function(reference_data, feature, grid, is_valid,
                                  predict_function, epsilon = 0.1,
                                  tolerance = 1e-7) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  validity <- validity_matrix(reference_data, feature, grid, is_valid)
  reference <- fit_relaxed_reference(validity, epsilon = epsilon, tolerance = tolerance)
  if (!reference$feasible) {
    stop(sprintf("relaxed reference constraints were not met (maximum violation %.3g).",
                 reference$max_violation), call. = FALSE)
  }
  predictions <- lapply(grid, function(z) {
    predictions <- predict_function(make_queries(reference_data, feature, z))
    if (!is.numeric(predictions) || length(predictions) != nrow(reference_data) ||
        any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per query row.", call. = FALSE)
    }
    predictions
  })
  prediction_matrix <- do.call(cbind, predictions)
  ordinary_pd <- colMeans(prediction_matrix)
  curve <- vapply(predictions, function(values) sum(reference$weights * values), numeric(1))
  density_ratio <- nrow(reference_data) * reference$weights
  reference_covariance <- colMeans(
    (density_ratio - 1) * sweep(prediction_matrix, 2L, ordinary_pd, "-")
  )
  rcpd_pd_difference <- curve - ordinary_pd
  weight_prediction_correlation <- vapply(predictions, function(values) {
    if (stats::sd(values) == 0 || stats::sd(reference$weights) == 0) 0
    else stats::cor(reference$weights, values)
  }, numeric(1))
  diagnostics <- reference$pattern_diagnostics
  data.frame(
    z = grid, relaxed_cspd = curve, ordinary_pd = ordinary_pd,
    rcpd_pd_difference = rcpd_pd_difference,
    reference_covariance = reference_covariance,
    centered_reference_covariance = reference_covariance - mean(reference_covariance),
    covariance_identity_error = abs(reference_covariance - rcpd_pd_difference),
    epsilon = epsilon,
    weighted_coverage = reference$coverage,
    min_coverage = min(reference$coverage), kl = reference$kl,
    mean_coverage = reference$mean_coverage,
    simultaneous_valid_mass = reference$simultaneous_valid_mass,
    effective_sample_size = reference$effective_sample_size,
    max_weight = max(reference$weights),
    weight_prediction_correlation = weight_prediction_correlation,
    observed_patterns = diagnostics$observed_patterns,
    median_pattern_size = diagnostics$median_pattern_size,
    min_pattern_size = diagnostics$min_pattern_size,
    singleton_pattern_fraction = diagnostics$singleton_pattern_fraction,
    active_constraints = diagnostics$active_constraints,
    active_affine_dimension = diagnostics$active_affine_dimension,
    active_licq = diagnostics$active_licq,
    binding_constraints = reference$binding_constraints,
    strict_complementarity = reference$strict_complementarity,
    epsilon_min = reference$epsilon_min,
    stringsAsFactors = FALSE
  )
}

#' Supremum shape difference after centering two curves on the same grid.
centered_curve_distance <- function(first, second) {
  if (!is.numeric(first) || !is.numeric(second) || length(first) != length(second) ||
      length(first) == 0L || any(!is.finite(first)) || any(!is.finite(second))) {
    stop("first and second must be finite numeric vectors of equal positive length.", call. = FALSE)
  }
  max(abs((first - mean(first)) - (second - mean(second))))
}

#' Centered sup-norm distance relative to the comparison curve's amplitude.
normalized_centered_curve_distance <- function(first, second,
                                               tolerance = sqrt(.Machine$double.eps)) {
  amplitude <- diff(range(second - mean(second)))
  if (!is.finite(amplitude) || amplitude <= tolerance) return(NA_real_)
  centered_curve_distance(first, second) / amplitude
}

#' Estimate common-support partial dependence (CSPD).
#'
#' @param predict_function A function mapping a data frame to one finite numeric
#'   prediction per row.
#' @return A data frame with CSPD, common-support mass and effective sample size.
estimate_cspd <- function(reference_data, feature, grid, is_valid, predict_function) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  common <- common_support_indicator(reference_data, feature, grid, is_valid)
  n_common <- sum(common)
  gamma <- mean(common)

  if (n_common == 0L) {
    return(data.frame(
      z = grid,
      cspd = NA_real_,
      gamma = gamma,
      n_common = n_common,
      stringsAsFactors = FALSE
    ))
  }

  cspd <- vapply(grid, function(z) {
    predictions <- predict_function(make_queries(reference_data[common, , drop = FALSE], feature, z))
    if (!is.numeric(predictions) || length(predictions) != n_common || any(!is.finite(predictions))) {
      stop("predict_function must return one finite numeric prediction per query row.", call. = FALSE)
    }
    mean(predictions)
  }, numeric(1))

  data.frame(z = grid, cspd = cspd, gamma = gamma, n_common = n_common, stringsAsFactors = FALSE)
}

#' Pointwise oracle Wald intervals for CSPD with a fixed validity rule.
#'
#' These intervals quantify only reference-sample variation conditional on a
#' fixed validity domain and fitted prediction function. They do not cover error
#' from estimating either object.
estimate_cspd_inference <- function(reference_data, feature, grid, is_valid, predict_function, conf_level = 0.95) {
  assert_data_frame(reference_data, "reference_data")
  assert_feature(reference_data, feature)
  grid <- assert_grid(grid)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must lie strictly between zero and one.", call. = FALSE)
  }
  common <- common_support_indicator(reference_data, feature, grid, is_valid)
  n_common <- sum(common)
  gamma <- mean(common)
  if (n_common == 0L) {
    return(data.frame(z = grid, cspd = NA_real_, std_error = NA_real_, lower = NA_real_, upper = NA_real_, gamma = gamma, n_common = n_common))
  }
  predictions <- lapply(grid, function(z) predict_function(make_queries(reference_data[common, , drop = FALSE], feature, z)))
  means <- vapply(predictions, mean, numeric(1))
  standard_errors <- vapply(predictions, function(value) {
    if (length(value) < 2L) NA_real_ else stats::sd(value) / sqrt(length(value))
  }, numeric(1))
  critical_value <- stats::qnorm(1 - (1 - conf_level) / 2)
  data.frame(
    z = grid, cspd = means, std_error = standard_errors,
    lower = means - critical_value * standard_errors,
    upper = means + critical_value * standard_errors,
    gamma = gamma, n_common = n_common, stringsAsFactors = FALSE
  )
}
