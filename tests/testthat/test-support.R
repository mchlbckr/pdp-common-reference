source(file.path("..", "..", "R", "support.R"))

testthat::test_that("ISC and its complement use the induced PDP queries", {
  reference_data <- data.frame(x = c(0, 1), u = c(0, 1))
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 0.1

  estimate <- estimate_isc(reference_data, "x", c(0, 1), valid)

  testthat::expect_equal(estimate$isc, c(0.5, 0.5))
  testthat::expect_equal(estimate$per, c(0.5, 0.5))
})

testthat::test_that("PDP averages every induced query over the reference rows", {
  reference_data <- data.frame(x = c(0, 1), u = c(2, 4))
  prediction <- function(newdata) newdata$x + newdata$u

  estimate <- estimate_pdp(reference_data, "x", c(0, 1), prediction)

  testthat::expect_equal(estimate$pdp, c(3, 4))
})

testthat::test_that("pointwise trimming exposes a changing reference population", {
  reference_data <- data.frame(x = c(0, 1, 2), u = c(0, 1, 2))
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 0.1
  prediction <- function(newdata) newdata$x + 10 * newdata$u

  estimate <- estimate_pointwise_trimmed_pdp(
    reference_data, "x", c(0, 2), valid, prediction
  )

  testthat::expect_equal(estimate$pointwise_trimmed_pdp, c(0, 22))
  testthat::expect_equal(estimate$n_valid, c(1L, 1L))
  testthat::expect_equal(estimate$valid_mass, c(1 / 3, 1 / 3))
})

testthat::test_that("range-only PD bounds have width equal to invalid mass times range", {
  reference_data <- data.frame(x = c(0, 1), u = c(0, 2))
  valid <- function(newdata) newdata$u == newdata$x
  prediction <- function(newdata) newdata$x + newdata$u
  bounds <- estimate_range_pd_bounds(
    reference_data, "x", 0, valid, prediction, c(-10, 10)
  )

  testthat::expect_equal(bounds$isc, 0.5)
  testthat::expect_equal(bounds$lower, -5)
  testthat::expect_equal(bounds$upper, 5)
  testthat::expect_equal(bounds$width, 20 * (1 - bounds$isc))
})

testthat::test_that("Lipschitz PD bounds sharpen range-only bounds and are sharp", {
  reference_data <- data.frame(x = c(0, 1), u = c(0, 2))
  valid <- function(newdata) newdata$u == newdata$x
  prediction <- function(newdata) newdata$x + newdata$u
  range_bounds <- estimate_range_pd_bounds(
    reference_data, "x", 0, valid, prediction, c(-10, 10)
  )
  lipschitz_bounds <- estimate_lipschitz_pd_bounds(
    reference_data, "x", 0, valid, prediction,
    lipschitz_constant = 1, prediction_bounds = c(-10, 10)
  )

  testthat::expect_equal(lipschitz_bounds$lower, -1)
  testthat::expect_equal(lipschitz_bounds$upper, 1)
  testthat::expect_lt(lipschitz_bounds$width, range_bounds$width)
  testthat::expect_true(lipschitz_bounds$lower <= 1 &
                          1 <= lipschitz_bounds$upper)
})

testthat::test_that("Lipschitz PD bounds reject an incompatible smoothness claim", {
  reference_data <- data.frame(x = c(0, 0), u = c(0, 1))
  valid <- function(newdata) rep(TRUE, nrow(newdata))
  prediction <- function(newdata) 2 * newdata$u

  testthat::expect_error(
    estimate_lipschitz_pd_bounds(
      reference_data, "x", 0, valid, prediction,
      lipschitz_constant = 1, prediction_bounds = c(-1, 3)
    ),
    "violate"
  )
})

testthat::test_that("CSPD reports when no common context exists", {
  reference_data <- data.frame(x = c(0, 1), u = c(0, 1))
  valid <- function(newdata) newdata$x == newdata$u
  prediction <- function(newdata) newdata$x + 10 * newdata$u

  estimate <- estimate_cspd(reference_data, "x", c(0, 1), valid, prediction)

  testthat::expect_equal(estimate$gamma, c(0, 0))
  testthat::expect_true(all(is.na(estimate$cspd)))
})

testthat::test_that("CSPD recovers additive contrasts on its common population", {
  reference_data <- data.frame(x = c(0, 1), u = c(0, 1))
  valid <- function(newdata) rep(TRUE, nrow(newdata))
  prediction <- function(newdata) newdata$x^2 + 10 * newdata$u

  estimate <- estimate_cspd(reference_data, "x", c(0, 1), valid, prediction)

  testthat::expect_equal(estimate$gamma, c(1, 1))
  testthat::expect_equal(estimate$cspd, c(5, 6))
  testthat::expect_equal(diff(estimate$cspd), 1)
})

testthat::test_that("ALE uses local changes and is centered", {
  reference_data <- data.frame(x = c(0.25, 0.75, 1.25, 1.75), u = c(2, 4, 6, 8))
  prediction <- function(newdata) newdata$x^2 + 100 * newdata$u

  estimate <- estimate_ale(reference_data, "x", c(0, 1, 2), prediction)

  testthat::expect_equal(estimate$n_bin, c(2L, 2L))
  # The curve is reported at interval midpoints, so its contrast is the
  # integral of 2x from 0.5 to 1.5, namely 2.
  testthat::expect_equal(diff(estimate$ale), 2)
  testthat::expect_equal(stats::weighted.mean(estimate$ale, estimate$n_bin), 0)
})

testthat::test_that("oracle CSPD intervals use the common reference rows", {
  reference_data <- data.frame(x = c(0, 1, 2), u = c(1, 2, 4))
  valid <- function(newdata) rep(TRUE, nrow(newdata))
  prediction <- function(newdata) newdata$x + newdata$u

  interval <- estimate_cspd_inference(reference_data, "x", c(0, 1), valid, prediction)

  testthat::expect_equal(interval$cspd, c(7 / 3, 10 / 3))
  testthat::expect_equal(interval$n_common, c(3L, 3L))
  testthat::expect_true(all(interval$lower < interval$cspd & interval$cspd < interval$upper))
})

testthat::test_that("conditional and subgroup baselines expose their populations", {
  reference_data <- data.frame(x = c(0, 0.1, 1, 1.1), u = c(1, 2, 3, 4))
  prediction <- function(newdata) newdata$x + newdata$u
  conditional <- estimate_conditional_curve(reference_data, "x", c(0, 1), prediction, bandwidth = 0.15)
  subgroup <- estimate_subgroup_pdp(reference_data, "x", c(0, 1), prediction,
                                    function(data) data$u <= 2)

  testthat::expect_equal(conditional$n_context, c(2, 2))
  testthat::expect_equal(subgroup$pdp, c(1.5, 2.5))
  testthat::expect_equal(subgroup$n_subgroup, c(2L, 2L))
  testthat::expect_equal(subgroup$subgroup_mass, c(0.5, 0.5))
})

testthat::test_that("relaxed KL projection attains simultaneous validity constraints", {
  validity <- matrix(c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE), nrow = 3)
  fit <- fit_relaxed_reference(validity, epsilon = 0.4)

  testthat::expect_true(fit$feasible)
  testthat::expect_true(all(fit$coverage >= 0.6 - 1e-7))
  testthat::expect_equal(sum(fit$weights), 1, tolerance = 1e-12)
  testthat::expect_gte(fit$effective_sample_size, 1)
  testthat::expect_gte(fit$kl, 0)
  testthat::expect_equal(fit$epsilon_min, 0, tolerance = 1e-7)
  testthat::expect_true(is.logical(fit$strict_complementarity))
  testthat::expect_true(fit$simultaneous_valid_mass >= 0 &&
                          fit$simultaneous_valid_mass <= 1)
  testthat::expect_equal(
    apply_relaxed_reference_tilt(validity, fit$multipliers),
    fit$weights, tolerance = 1e-12
  )
  boundary_multipliers <- fit$multipliers
  boundary_multipliers[[which.min(boundary_multipliers)]] <- -1e-18
  testthat::expect_equal(
    apply_relaxed_reference_tilt(validity, boundary_multipliers),
    apply_relaxed_reference_tilt(validity, pmax(boundary_multipliers, 0)),
    tolerance = 1e-12
  )
})

testthat::test_that("single-constraint references distinguish average and intersection validity", {
  validity <- rbind(
    c(FALSE, TRUE, TRUE),
    c(TRUE, FALSE, TRUE),
    c(TRUE, TRUE, FALSE)
  )
  average <- fit_average_validity_reference(validity, epsilon = 1 / 3)

  testthat::expect_true(average$feasible)
  testthat::expect_equal(average$constraint_value, 2 / 3, tolerance = 1e-8)
  testthat::expect_equal(average$simultaneous_valid_mass, 0)
  testthat::expect_error(
    fit_intersection_reference(validity, epsilon = 0.2),
    "epsilon >= 1",
    fixed = TRUE
  )
})

testthat::test_that("intersection reference directly controls common-valid mass", {
  validity <- rbind(
    c(TRUE, TRUE), c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE)
  )
  fit <- fit_intersection_reference(validity, epsilon = 0.25)

  testthat::expect_true(fit$feasible)
  testthat::expect_equal(fit$simultaneous_valid_mass, 0.75, tolerance = 1e-7)
  testthat::expect_equal(fit$constraint_value, 0.75, tolerance = 1e-7)
  testthat::expect_true(all(fit$coverage >= 0.75 - 1e-7))
})

testthat::test_that("weighted PDP uses one fixed reference distribution", {
  reference_data <- data.frame(x = c(0, 1), u = c(2, 4))
  prediction <- function(newdata) newdata$x + newdata$u
  estimate <- estimate_weighted_pdp(
    reference_data, "x", c(0, 1), prediction, c(0.25, 0.75),
    value_name = "curve"
  )

  testthat::expect_equal(estimate$curve, c(3.5, 4.5))
})

testthat::test_that("adaptive constraint generation controls a dense grid", {
  reference_data <- data.frame(
    x = rep(0, 101), u = seq(-1, 1, length.out = 101)
  )
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 1
  fit <- fit_adaptive_relaxed_reference(
    reference_data, "x", c(-0.8, 0.8), valid,
    epsilon = 0.25, initial_grid = c(-0.8, 0, 0.8),
    validation_grid = seq(-0.8, 0.8, length.out = 41)
  )

  testthat::expect_true(fit$converged)
  testthat::expect_gte(fit$min_validation_coverage, 0.75 - 2e-5)
  testthat::expect_true(all(fit$active_grid %in% fit$validation_grid))
  testthat::expect_equal(sum(fit$weights), 1, tolerance = 1e-10)
})

testthat::test_that("critical interval cells attain every weighted coverage minimum", {
  lower <- c(-Inf, -0.4, 0.1, 0.6)
  upper <- c(-0.2, 0.5, 0.8, Inf)
  interval <- c(-1, 1)
  critical <- interval_critical_grid(lower, upper, interval)
  critical_patterns <- interval_validity_matrix(lower, upper, critical)
  audit <- seq(interval[[1]], interval[[2]], length.out = 2001L)
  audit_patterns <- interval_validity_matrix(lower, upper, audit)

  weight_sets <- list(
    rep(1 / length(lower), length(lower)),
    c(0.1, 0.2, 0.3, 0.4),
    c(0.55, 0.15, 0.15, 0.15)
  )
  for (weights in weight_sets) {
    testthat::expect_equal(
      min(as.numeric(crossprod(weights, critical_patterns))),
      min(as.numeric(crossprod(weights, audit_patterns))),
      tolerance = 1e-12
    )
  }
  testthat::expect_false(any(c(interval, -0.4, -0.2, 0.1, 0.5, 0.6, 0.8) %in%
                               critical))
})

testthat::test_that("interval exchange solver equals the exact continuum reduction", {
  lower <- c(-0.9, -0.7, -0.2, 0.1, 0.4, 0.65)
  upper <- c(-0.1, 0.2, 0.55, 0.8, 0.95, 1.1)
  interval <- c(-0.75, 0.75)
  epsilon <- 0.7
  exact <- fit_interval_relaxed_reference(
    lower, upper, interval, epsilon = epsilon, tolerance = 2e-5
  )
  adaptive <- fit_adaptive_interval_relaxed_reference(
    lower, upper, interval, epsilon = epsilon,
    initial_grid = c(-0.75, 0, 0.75), tolerance = 2e-5
  )

  testthat::expect_true(adaptive$converged)
  testthat::expect_gte(adaptive$min_critical_coverage,
                       1 - epsilon - 2e-5)
  testthat::expect_equal(adaptive$weights, exact$weights, tolerance = 5e-5)
  testthat::expect_equal(adaptive$kl, exact$kl, tolerance = 5e-5)
  testthat::expect_lte(adaptive$iterations, length(exact$critical_grid))
})

testthat::test_that("relaxed-reference range bounds are controlled by epsilon", {
  reference_data <- data.frame(x = 0, u = c(-1, 0, 1, 2))
  grid <- c(0, 1)
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 1
  prediction <- function(newdata) pmin(3, pmax(-2, newdata$x + newdata$u))
  validity <- validity_matrix(reference_data, "x", grid, valid)
  epsilon <- 0.5
  fit <- fit_relaxed_reference(validity, epsilon = epsilon)
  bounds <- estimate_range_pd_bounds(
    reference_data, "x", grid, valid, prediction,
    prediction_bounds = c(-2, 3), weights = fit$weights
  )

  testthat::expect_equal(bounds$isc, fit$coverage, tolerance = 1e-8)
  testthat::expect_lte(max(bounds$width), 5 * epsilon + 1e-7)
})

testthat::test_that("conditional curve bootstrap returns simultaneous bands", {
  reference_data <- data.frame(x = rep(0, 40), u = seq(-1, 1, length.out = 40))
  grid <- c(-0.5, 0, 0.5)
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 1
  prediction <- function(newdata) newdata$x^2 + newdata$u
  h <- validity_matrix(reference_data, "x", grid, valid)
  p <- prediction_matrix(reference_data, "x", grid, prediction)
  fit <- bootstrap_reference_curve_stability(
    h, p, grid, epsilon = 0.25, repetitions = 10L, seed = 1,
    workers = 1L
  )

  testthat::expect_equal(sort(unique(fit$bands$method)),
                         c("Hard CSPD", "PDP", "RCPD"))
  testthat::expect_equal(nrow(fit$metrics), 10L)
  testthat::expect_true(all(fit$bands$lower <= fit$bands$estimate &
                              fit$bands$estimate <= fit$bands$upper))
})

testthat::test_that("minimum epsilon detects hard feasibility and incompatible patterns", {
  hard_feasible <- rbind(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE))
  disjoint <- rbind(c(TRUE, FALSE), c(FALSE, TRUE))

  testthat::expect_equal(minimum_relaxation_epsilon(hard_feasible), 0, tolerance = 1e-8)
  testthat::expect_equal(minimum_relaxation_epsilon(disjoint), 0.5, tolerance = 1e-8)
  testthat::expect_error(
    fit_relaxed_reference(disjoint, epsilon = 0.4),
    "epsilon >= 0.5",
    fixed = TRUE
  )
})

testthat::test_that("validity pattern diagnostics report cells and active constraints", {
  validity <- rbind(
    c(TRUE, TRUE), c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE)
  )
  diagnostics <- validity_pattern_diagnostics(
    validity, weights = c(0.4, 0.3, 0.2, 0.1), multipliers = c(0.2, 0)
  )

  testthat::expect_equal(diagnostics$observed_patterns, 3)
  testthat::expect_equal(diagnostics$median_pattern_size, 1)
  testthat::expect_equal(diagnostics$min_pattern_size, 1)
  testthat::expect_equal(diagnostics$singleton_pattern_fraction, 2 / 3)
  testthat::expect_equal(diagnostics$active_constraints, 1)
  testthat::expect_equal(diagnostics$active_affine_dimension, 1)
  testthat::expect_true(diagnostics$active_licq)
  testthat::expect_equal(diagnostics$max_weight, 0.4)
  testthat::expect_equal(diagnostics$epsilon_min, 0)

  rank_deficient <- validity_pattern_diagnostics(
    rbind(c(TRUE, TRUE), c(FALSE, FALSE)),
    multipliers = c(0.2, 0.3)
  )
  testthat::expect_equal(rank_deficient$active_affine_dimension, 1)
  testthat::expect_false(rank_deficient$active_licq)
})

testthat::test_that("inactive relaxed constraints preserve the empirical reference", {
  validity <- matrix(c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), nrow = 3)
  fit <- fit_relaxed_reference(validity, epsilon = 0.7)

  testthat::expect_true(fit$feasible)
  testthat::expect_equal(fit$weights, rep(1 / 3, 3), tolerance = 1e-8)
  testthat::expect_equal(fit$kl, 0, tolerance = 1e-8)
})

testthat::test_that("epsilon zero recovers hard CSPD exactly", {
  reference_data <- data.frame(x = c(0, 1, 2), u = c(0, 1, 2))
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 1
  prediction <- function(newdata) newdata$x + 2 * newdata$u
  grid <- c(0, 1)
  hard <- estimate_cspd(reference_data, "x", grid, valid, prediction)
  relaxed <- estimate_relaxed_cspd(reference_data, "x", grid, valid, prediction, epsilon = 0)

  testthat::expect_equal(relaxed$relaxed_cspd, hard$cspd)
  testthat::expect_equal(relaxed$weighted_coverage, rep(1, length(grid)))
  testthat::expect_equal(min(relaxed$weighted_coverage), unique(relaxed$min_coverage))
  testthat::expect_equal(relaxed$effective_sample_size, hard$n_common)
  testthat::expect_equal(relaxed$kl, -log(hard$gamma))
  testthat::expect_equal(unique(relaxed$epsilon_min), 0)
})

testthat::test_that("relaxed-reference covariance identity holds gridwise", {
  reference_data <- data.frame(x = c(0, 1, 2, 3), u = c(-1, 0, 1, 2))
  grid <- c(0, 1)
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 1
  prediction <- function(newdata) newdata$x + newdata$x * newdata$u + newdata$u^2
  relaxed <- estimate_relaxed_cspd(
    reference_data, "x", grid, valid, prediction, epsilon = 0.5
  )

  testthat::expect_equal(
    relaxed$reference_covariance,
    relaxed$rcpd_pd_difference,
    tolerance = 1e-10
  )
  testthat::expect_lte(max(relaxed$covariance_identity_error), 1e-10)
  testthat::expect_equal(
    relaxed$centered_reference_covariance,
    relaxed$reference_covariance - mean(relaxed$reference_covariance),
    tolerance = 1e-12
  )
})

testthat::test_that("centered curve distance removes pure level shifts", {
  testthat::expect_equal(centered_curve_distance(c(1, 2, 4), c(11, 12, 14)), 0)
  testthat::expect_gt(centered_curve_distance(c(1, 2, 4), c(1, 3, 4)), 0)
  testthat::expect_equal(
    normalized_centered_curve_distance(c(1, 3, 4), c(1, 2, 4)),
    2 / 9
  )
})

testthat::test_that("common-reference estimators preserve additive curve shape", {
  reference_data <- data.frame(x = c(-1, 0, 1, 2), u = c(-1, 0, 1, 3))
  grid <- c(-0.5, 0.5, 1.5)
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 2
  prediction <- function(newdata) newdata$x^2 + 3 * newdata$u

  pdp <- estimate_pdp(reference_data, "x", grid, prediction)$pdp
  hard <- estimate_cspd(reference_data, "x", grid, valid, prediction)$cspd
  relaxed <- estimate_relaxed_cspd(
    reference_data, "x", grid, valid, prediction, epsilon = 0.4
  )$relaxed_cspd

  testthat::expect_equal(centered_curve_distance(pdp, hard), 0, tolerance = 1e-12)
  testthat::expect_equal(centered_curve_distance(pdp, relaxed), 0, tolerance = 1e-12)
})

testthat::test_that("RCPD reports pattern and weight-prediction diagnostics", {
  reference_data <- data.frame(x = c(0, 1, 2, 3), u = c(0, 1, 2, 3))
  grid <- c(1, 2)
  valid <- function(newdata) abs(newdata$x - newdata$u) <= 2
  prediction <- function(newdata) newdata$x * newdata$u
  estimate <- estimate_relaxed_cspd(
    reference_data, "x", grid, valid, prediction, epsilon = 0.25
  )

  testthat::expect_true(all(is.finite(estimate$weight_prediction_correlation)))
  testthat::expect_gte(unique(estimate$observed_patterns), 1)
  testthat::expect_gte(unique(estimate$active_constraints), 0)
  testthat::expect_gte(unique(estimate$binding_constraints), 0)
  testthat::expect_equal(unique(estimate$epsilon_min), 0)
})
