source(file.path("..", "..", "R", "support.R"))
source(file.path("..", "..", "simulations", "s1_additive_correlation.R"))
source(file.path("..", "..", "simulations", "s4_off_support_perturbation.R"))
source(file.path("..", "..", "R", "validity.R"))
source(file.path("..", "..", "R", "plot_theme.R"))
source(file.path("..", "..", "simulations", "s7_mtcars_case_study.R"))
source(file.path("..", "..", "simulations", "s2_nonlinear_manifold.R"))
source(file.path("..", "..", "simulations", "s8_oracle_inference.R"))
source(file.path("..", "..", "simulations", "s9_estimated_validity_inference.R"))
source(file.path("..", "..", "simulations", "s10_outer_bootstrap.R"))
source(file.path("..", "..", "simulations", "s11_bootstrap_sample_size.R"))
source(file.path("..", "..", "simulations", "s12_validity_learner_robustness.R"))
source(file.path("..", "..", "simulations", "s13_airquality_case_study.R"))
source(file.path("..", "..", "simulations", "s14_experiment_matrix.R"))
source(file.path("..", "..", "simulations", "s16_wine_quality_case_study.R"))
source(file.path("..", "..", "simulations", "s17_large_mixed_case_studies.R"))
source(file.path("..", "..", "simulations", "s15_boston_case_study.R"))
source(file.path("..", "..", "simulations", "s18_relaxed_reference.R"))
source(file.path("..", "..", "simulations", "s19_calibration_convergence.R"))
source(file.path("..", "..", "simulations", "s20_alpha_sensitivity.R"))
source(file.path("..", "..", "simulations", "s21_dimension_sensitivity.R"))
source(file.path("..", "..", "simulations", "s22_rcpd_rate.R"))
source(file.path("..", "..", "simulations", "s23_composition_trap.R"))
source(file.path("..", "..", "simulations", "s24_constraint_comparison.R"))
source(file.path("..", "..", "simulations", "s25_error_decomposition.R"))
source(file.path("..", "..", "simulations", "s26_grid_resolution.R"))
source(file.path("..", "..", "simulations", "s27_application_uncertainty.R"))
source(file.path("..", "..", "simulations", "s28_full_split_stability.R"))
source(file.path("..", "..", "simulations", "s29_rcpd_comparison.R"))
source(file.path("..", "..", "simulations", "s30_lipschitz_wine.R"))
source(file.path("..", "..", "simulations", "s31_projection_generalization.R"))
source(file.path("..", "..", "simulations", "s32_bank_offsupport_stress.R"))
source(file.path("..", "..", "simulations", "s34_conditional_subgroup_comparison.R"))
local({
  previous_directory <- getwd()
  on.exit(setwd(previous_directory), add = TRUE)
  setwd(file.path("..", ".."))
  source(file.path("simulations", "s33_california_longitude_application.R"),
         local = globalenv())
})

testthat::test_that("S4 models agree exactly on valid queries", {
  data <- simulate_s4_data(n = 200L, rho = 0.7, seed = 1)
  valid <- s1_validity_rule(rho = 0.7)
  query_data <- make_queries(data, "x1", 0.25)
  on_validity_domain <- valid(query_data)

  testthat::expect_equal(
    s4_base_model(query_data)[on_validity_domain],
    s4_perturbed_model(query_data, rho = 0.7)[on_validity_domain],
    tolerance = 1e-14
  )
})

testthat::test_that("S4 respects the empirical support sensitivity bound", {
  simulation <- run_s4(n = 1000L, seed = 2)

  testthat::expect_true(all(simulation$pdp_gap <= simulation$bound + 1e-12))
  testthat::expect_gt(unique(simulation$gamma), 0)
  testthat::expect_lt(max(simulation$cspd_gap), 1e-12)
})

testthat::test_that("S7 reports feasible common support and grid sensitivity", {
  case_study <- run_s7(grid_points = 21L)

  testthat::expect_gt(case_study$summary$gamma, 0)
  testthat::expect_lt(case_study$summary$gamma, 1)
  testthat::expect_equal(nrow(case_study$curves), 21L)
  # Narrower grids relax the common-context requirement, so their mass cannot
  # decrease as the central-quantile width shrinks.
  testthat::expect_true(all(diff(case_study$grid_sensitivity$gamma) >= -1e-12))
})

testthat::test_that("S8 oracle intervals have reasonable finite-sample coverage", {
  simulation <- run_s8(n = 250L, repetitions = 100L, grid = seq(-0.5, 0.5, length.out = 7L), seed = 4)

  testthat::expect_true(all(simulation$coverage > 0.75))
  testthat::expect_true(all(simulation$coverage < 1.01))
})

testthat::test_that("S9 separates plug-in and oracle targets after sample splitting", {
  simulation <- run_s9(repetitions = 80L, grid = seq(-0.5, 0.5, length.out = 7L), seed = 5)

  testthat::expect_true(all(simulation$covered_plugin > 0.7))
  testthat::expect_true(all(simulation$gamma_plugin > 0))
  testthat::expect_true(all(simulation$gamma_evaluation > 0))
})

testthat::test_that("S10 outer bootstrap returns finite oracle-target intervals", {
  simulation <- run_s10(repetitions = 20L, bootstrap_repetitions = 30L,
                        grid = seq(-0.5, 0.5, length.out = 5L), seed = 6)

  testthat::expect_true(all(is.finite(simulation$bootstrap_width)))
  testthat::expect_true(all(simulation$covered_bootstrap >= 0 & simulation$covered_bootstrap <= 1))
})

testthat::test_that("S11 summarizes bootstrap calibration by sample size", {
  simulation <- run_s11(split_sizes = c(60L, 100L), repetitions = 10L,
                        bootstrap_repetitions = 10L, grid = seq(-0.5, 0.5, length.out = 5L), seed = 7)

  testthat::expect_equal(nrow(simulation), 2L)
  testthat::expect_true(all(is.finite(simulation$bootstrap_width)))
})

testthat::test_that("S12 compares validity learners on nonlinear support", {
  simulation <- run_s12(n_training = 80L, n_calibration = 80L, n_evaluation = 150L,
                        repetitions = 4L, k = 8L, grid = seq(-0.3, 0.3, length.out = 7L), seed = 8)

  testthat::expect_equal(nrow(simulation), 2L)
  testthat::expect_true(all(is.finite(simulation$isc_mae)))
  testthat::expect_true(all(simulation$gamma >= 0 & simulation$gamma <= 1))
})

testthat::test_that("S13 supplies a nonempty observational common population", {
  case_study <- run_s13(grid_points = 21L)

  testthat::expect_gt(unique(case_study$gamma), 0)
  testthat::expect_lt(unique(case_study$gamma), 1)
  testthat::expect_equal(nrow(case_study), 21L)
})

testthat::test_that("S14 runs baselines and all validity learners in one protocol", {
  simulation <- run_s14(rhos = 0.5, grid_half_widths = 0.4, model_classes = "linear",
                        repetitions = 1L, n_training = 80L, n_calibration = 80L,
                        n_evaluation = 120L, k = 10L, seed = 9)
  summary <- summarize_s14(simulation)

  testthat::expect_true(all(c("PDP", "ALE", "conditional curve", "CSPD") %in% simulation$method))
  testthat::expect_equal(sort(na.omit(unique(simulation$learner))),
                         c("conditional_density", "knn", "locally_scaled", "oracle", "split_conformal"))
  testthat::expect_true(all(is.finite(simulation$shape_distance_to_oracle)))
  testthat::expect_true(all(summary$n_repetitions == 1L))
  testthat::expect_true(all(simulation$residual_density_identical))
})

testthat::test_that("S14 residual and homoskedastic density rules induce identical H matrices", {
  set.seed(901)
  training <- simulate_s1_data(120L, rho = 0.6)
  calibration <- simulate_s1_data(100L, rho = 0.6)
  evaluation <- simulate_s1_data(90L, rho = 0.6)
  grid <- seq(-0.8, 0.8, length.out = 17L)
  residual_fit <- fit_split_conformal_validity(training, calibration, "x1", alpha = 0.1)
  density_fit <- fit_conditional_density_validity(training, calibration, "x1", alpha = 0.1)

  residual_h <- validity_matrix(evaluation, "x1", grid, function(newdata) {
    predict_validity(residual_fit, newdata)
  })
  density_h <- validity_matrix(evaluation, "x1", grid, function(newdata) {
    predict_conditional_density_validity(density_fit, newdata)
  })

  testthat::expect_identical(residual_h, density_h)
})

testthat::test_that("S15 uses a substantial housing benchmark without the race-derived field", {
  case_study <- run_s15(grid_points = 21L)

  testthat::expect_equal(unique(case_study$n_common), as.integer(unique(case_study$gamma) * 506))
  testthat::expect_gt(unique(case_study$gamma), 0)
  testthat::expect_lt(unique(case_study$gamma), 1)
})

testthat::test_that("S16 keeps model, validity, and evaluation data separate", {
  case_study <- run_s16(grid_points = 21L, seed = 10)

  testthat::expect_equal(nrow(case_study), 42L)
  testthat::expect_equal(length(unique(case_study$model_case)), 2L)
  testthat::expect_true("additive_gam" %in% case_study$model_name)
  testthat::expect_true(attr(case_study, "selected_nonadditive_model") %in%
                          c("interaction_gam", "boosted_tree"))
  testthat::expect_equal(nrow(attr(case_study, "candidate_summary")), 3L)
  testthat::expect_equal(nrow(attr(case_study, "evaluation_contexts")),
                         unique(case_study$n_evaluation))
  testthat::expect_equal(sort(unique(attr(case_study, "support_profile")$alpha)),
                         c(0.01, 0.05, 0.10, 0.20))
  testthat::expect_gt(unique(case_study$n_common), 0)
  testthat::expect_true(all(is.finite(case_study$heldout_mae)))
  testthat::expect_equal(unique(case_study$n_training + case_study$n_calibration + case_study$n_evaluation), 1599)
  testthat::expect_gte(unique(case_study$observed_patterns), 1)
  testthat::expect_lte(unique(case_study$epsilon_min), unique(case_study$relaxed_epsilon))
})

testthat::test_that("S17 shopper application retains an independent common population", {
  case_study <- run_s17_shoppers(seed = 11)

  testthat::expect_gt(unique(case_study$n_common), 0)
  testthat::expect_true(all(case_study$isc >= 0 & case_study$isc <= 1))
  testthat::expect_true(all(is.finite(case_study$heldout_brier)))
  testthat::expect_true(all(is.finite(case_study$weight_prediction_correlation)))
  testthat::expect_gte(unique(case_study$active_constraints), 0)
  testthat::expect_equal(sort(unique(attr(case_study, "support_profile")$alpha)),
                         c(0.01, 0.05, 0.10, 0.20))
})

testthat::test_that("S17 bank application compares conditional and joint scores", {
  case_study <- run_s17_bank(seed = 12)
  comparison <- attr(case_study, "score_comparison")

  testthat::expect_equal(nrow(comparison), 2L)
  testthat::expect_true(all(comparison$mean_isc >= 0 & comparison$mean_isc <= 1))
  testthat::expect_true(all(comparison$gamma >= 0 & comparison$gamma <= 1))
  testthat::expect_equal(unique(case_study$joint_gamma), comparison$gamma[[2]])
})

testthat::test_that("S18 exposes the hard-to-relaxed reference trade-off", {
  simulation <- run_s18(repetitions = 3L, n_evaluation = 300L,
                        epsilons = c(0, 0.1, 0.2), seed = 12, workers = 1L)
  summary <- summarize_s18(simulation)

  testthat::expect_equal(nrow(summary), 3L)
  testthat::expect_true(all(summary$min_coverage >= 1 - summary$epsilon - 2e-5))
  testthat::expect_true(all(diff(summary$kl) <= 1e-5))
  testthat::expect_true(all(diff(summary$effective_sample_fraction) >= -1e-5))
  testthat::expect_true(all(summary$epsilon_min <= summary$epsilon + 2e-5))
  testthat::expect_true(all(summary$observed_patterns >= 1))
})

testthat::test_that("S19 returns calibration convergence diagnostics", {
  simulation <- run_s19(split_sizes = c(60L, 120L), repetitions = 3L,
                        n_evaluation = 300L, seed = 13, workers = 1L)
  summary <- summarize_s19(simulation)

  testthat::expect_equal(nrow(summary), 2L)
  testthat::expect_true(all(summary$learned_gamma >= 0 & summary$learned_gamma <= 1))
  testthat::expect_true(all(summary$radius_ratio > 0))
})

testthat::test_that("S20 reports alpha sensitivity on a common oracle scale", {
  simulation <- run_s20(alphas = c(0.05, 0.20), repetitions = 3L,
                        n_evaluation = 300L, seed = 14, workers = 1L)
  summary <- summarize_s20(simulation)

  testthat::expect_equal(nrow(summary), 2L)
  testthat::expect_true(all(summary$min_isc <= summary$mean_isc))
  testthat::expect_true(all(summary$mean_isc <= summary$max_isc))
  testthat::expect_gt(summary$gamma[[1]], summary$gamma[[2]])
})

testthat::test_that("S21 generates the requested one-factor context dimensions", {
  generated <- simulate_s21_data(30L, context_dimension = 4L, rho = 0.5, seed = 15)

  testthat::expect_equal(dim(generated), c(30L, 5L))
  testthat::expect_equal(names(generated), c("x1", "u1", "u2", "u3", "u4"))
})

testthat::test_that("S21 oracle validity functions return one logical per query", {
  generated <- simulate_s21_data(25L, context_dimension = 4L, rho = 0.5, seed = 16)
  conditional <- s21_oracle_validity(4L, 0.5, 0.05, "conditional")
  joint <- s21_oracle_validity(4L, 0.5, 0.05, "joint")

  testthat::expect_type(conditional(generated), "logical")
  testthat::expect_type(joint(generated), "logical")
  testthat::expect_length(conditional(generated), 25L)
  testthat::expect_length(joint(generated), 25L)
})

testthat::test_that("S21 runs, summarizes, and plots dimension sensitivity", {
  simulation <- run_s21(
    context_dimensions = c(1L, 4L), repetitions = 3L, n_evaluation = 120L,
    grid = seq(-0.8, 0.8, length.out = 9L), seed = 17, workers = 1L
  )
  summary <- summarize_s21(simulation)
  output <- tempfile(fileext = ".png")
  plotted <- plot_s21(summary, output)

  testthat::expect_equal(nrow(simulation), 12L)
  testthat::expect_equal(nrow(summary), 4L)
  testthat::expect_true(all(summary$epsilon_min <= summary$fitted_epsilon + 1e-8))
  testthat::expect_length(plotted, 4L)
  testthat::expect_true(all(file.exists(plotted)))
})

testthat::test_that("S22 builds a finite population-approximation target", {
  target <- s22_oracle_target(
    n_oracle = 1000L, grid = seq(-0.6, 0.6, length.out = 7L), seed = 18
  )

  testthat::expect_length(target, 7L)
  testthat::expect_true(all(is.finite(target)))
})

testthat::test_that("S22 runs, summarizes, and plots the RCPD rate diagnostic", {
  simulation <- run_s22(
    sample_sizes = c(120L, 240L, 480L), repetitions = 3L,
    grid = seq(-0.6, 0.6, length.out = 7L), n_oracle = 1500L,
    seed = 19, workers = 1L
  )
  summary <- summarize_s22(simulation)
  output <- tempfile(fileext = ".png")
  plotted <- plot_s22(summary, output)

  testthat::expect_equal(nrow(simulation), 9L)
  testthat::expect_equal(nrow(summary), 3L)
  testthat::expect_true(is.finite(unique(summary$log_log_slope)))
  testthat::expect_true(file.exists(plotted))
})

testthat::test_that("S23 isolates pointwise composition bias", {
  simulation <- run_s23(
    repetitions = 3L, n_evaluation = 300L,
    grid = seq(-0.5, 0.5, length.out = 9L), seed = 20,
    workers = 1L
  )
  summary <- summarize_s23(simulation)
  hard <- summary[summary$method == "Hard CSPD", ]
  trimmed <- summary[summary$method == "Pointwise trim", ]
  output <- tempfile(fileext = ".png")

  testthat::expect_equal(nrow(simulation$metrics), 3L * 2L * 5L)
  testthat::expect_length(simulation$sample_focal, 300L)
  testthat::expect_true(all(hard$normalized_shape_error <
                              trimmed$normalized_shape_error))
  testthat::expect_true(all(summary$gamma > 0))
  plotted <- plot_s23(simulation, output)
  testthat::expect_length(plotted, 4L)
  testthat::expect_true(all(file.exists(plotted)))
})

testthat::test_that("S24 distinguishes pointwise, average, and intersection constraints", {
  simulation <- run_s24(
    repetitions = 3L, n_evaluation = 400L, k = 6L,
    epsilon = 0.18, seed = 21, workers = 1L
  )
  summary <- summarize_s24(simulation)
  pointwise <- summary$diagnostics[
    summary$diagnostics$reference == "Pointwise constraints", ]
  average <- summary$diagnostics[
    summary$diagnostics$reference == "Average constraint", ]
  intersection <- summary$diagnostics[
    summary$diagnostics$reference == "Intersection constraint", ]
  output_directory <- tempfile()

  testthat::expect_equal(unique(simulation$metrics$gamma), 0)
  testthat::expect_true(pointwise$min_coverage >= 0.82 - 2e-5)
  testthat::expect_lt(average$min_coverage, pointwise$min_coverage)
  testthat::expect_equal(intersection$feasible_rate, 0)
  output <- plot_s24(simulation, output_directory)
  testthat::expect_length(output, 4L)
  testthat::expect_true(all(file.exists(output)))
})

testthat::test_that("S25 separates evaluation and validity-target error", {
  simulation <- run_s25(
    split_sizes = c(80L, 160L), repetitions = 5L,
    grid = seq(-0.5, 0.5, length.out = 5L), seed = 22,
    workers = 1L
  )
  summary <- summarize_s25(simulation)
  output <- tempfile(fileext = ".png")

  testthat::expect_equal(nrow(summary), 2L)
  testthat::expect_true(all(summary$evaluation_rmse > 0))
  testthat::expect_true(all(summary$validity_rmse > 0))
  testthat::expect_true(all(summary$plugin_coverage >= 0 &
                              summary$plugin_coverage <= 1))
  testthat::expect_true(all(file.exists(plot_s25(summary, output))))
})

testthat::test_that("S26 verifies endpoint identity and exact continuum coverage", {
  simulation <- run_s26(
    repetitions = 3L, n_evaluation = 300L,
    grid_sizes = c(3L, 7L, 15L), validation_size = 31L,
    seed = 23, workers = 1L
  )
  summary <- summarize_s26(simulation)
  fixed <- simulation[simulation$method == "Fixed grid", ]
  exact <- simulation[simulation$method == "Exact exchange", ]
  output <- tempfile(fileext = ".png")

  by_rep <- split(fixed$hard_gamma, fixed$repetition)
  testthat::expect_true(all(vapply(by_rep, function(x) length(unique(x)) == 1L,
                                    logical(1))))
  testthat::expect_true(all(exact$min_continuum_coverage >= 0.9 - 2e-5))
  testthat::expect_true(all(exact$curve_distance_to_exact == 0))
  testthat::expect_true(all(exact$constraints_used <=
                              2L * 300L + 3L))
  testthat::expect_true(all(file.exists(plot_s26(summary, output))))
})

testthat::test_that("S27 reports conditional application stability bands", {
  wine <- run_s16(grid_points = 9L, seed = 24)
  shoppers <- run_s17_shoppers(seed = 25)
  simulation <- run_s27(
    wine, list(shoppers), repetitions = 5L, seed = 26, workers = 1L
  )
  summary <- summarize_s27(simulation)
  output <- tempfile(fileext = ".png")

  testthat::expect_equal(nrow(summary), 3L)
  testthat::expect_true(all(simulation$differences$lower <=
                              simulation$differences$estimate))
  testthat::expect_true(all(simulation$differences$estimate <=
                              simulation$differences$upper))
  testthat::expect_true(all(file.exists(plot_s27(simulation, output))))
})

testthat::test_that("S28 summarizes full-pipeline split stability", {
  simulation <- run_s28(
    wine_repetitions = 2L, bank_repetitions = 1L,
    shopper_repetitions = 1L, seed = 27, workers = 1L
  )
  summary <- summarize_s28(simulation)
  output <- tempfile(fileext = ".png")

  testthat::expect_equal(nrow(simulation), 6L)
  testthat::expect_equal(nrow(summary), 4L)
  testthat::expect_true(all(summary$median_normalized_distance >= 0))
  testthat::expect_true(all(file.exists(plot_s28(simulation, output))))
})

testthat::test_that("S29 compares all methods against one external target", {
  simulation <- run_s29(
    repetitions = 12L, n_evaluation = 400L, seed = 28, workers = 1L
  )
  summary <- summarize_s29(simulation)
  errors <- stats::setNames(
    summary$metrics$normalized_shape_error, summary$metrics$method
  )
  output <- tempfile(fileext = ".png")

  testthat::expect_equal(
    sort(summary$metrics$method),
    sort(c("Ordinary PD", "PTPD", "Hard CSPD", "RCPD", "ALE"))
  )
  testthat::expect_lt(errors[["RCPD"]], errors[["Hard CSPD"]])
  testthat::expect_lt(errors[["RCPD"]], errors[["PTPD"]])
  testthat::expect_lt(errors[["RCPD"]], errors[["ALE"]])
  testthat::expect_true(all(file.exists(plot_s29(summary, output))))
})

testthat::test_that("S30 produces nested admissible Lipschitz sensitivity bands", {
  wine <- run_s16(grid_points = 9L, seed = 29)
  simulation <- run_s30(wine, multipliers = c(1, 1.5, 2))
  summary <- summarize_s30(simulation)
  finite_widths <- summary$widths[
    summary$widths$scenario != "Range only", "mean_width"
  ]
  output <- tempfile(fileext = ".png")

  testthat::expect_gt(summary$admissible_minimum, 0)
  testthat::expect_equal(summary$model_name, "interaction_gam")
  testthat::expect_true(is.logical(summary$pdp_inside_lf0))
  testthat::expect_true(is.finite(summary$minimum_pdp_band_margin))
  testthat::expect_true(is.finite(summary$maximum_pdp_band_violation))
  testthat::expect_gte(summary$maximum_pdp_band_violation, 0)
  testthat::expect_true(all(diff(finite_widths) >= -1e-10))
  testthat::expect_lt(max(finite_widths),
                      summary$widths$mean_width[
                        summary$widths$scenario == "Range only"
                      ])
  testthat::expect_true(all(file.exists(plot_s30(simulation, output))))
})

testthat::test_that("S31 certifies fitted tilts on independent contexts", {
  simulation <- run_s31(
    repetitions = 4L, audit_sizes = c(500L, 1000L),
    bounds = c("empirical_bernstein", "hoeffding"),
    n_projection = 300L, fitted_epsilons = c(0.02, 0.10),
    grid = seq(-0.8, 0.8, length.out = 9L), truth_size = 5000L,
    seed = 30, workers = 1L
  )
  summary <- summarize_s31(simulation)

  testthat::expect_equal(nrow(summary), 4L)
  testthat::expect_true(all(simulation$naive_truth_min_coverage >= 0 &
                              simulation$naive_truth_min_coverage <= 1))
  testthat::expect_true(all(summary$certification_rate >= 0 &
                              summary$certification_rate <= 1))
  selected <- simulation[simulation$certified, , drop = FALSE]
  if (nrow(selected)) {
    testthat::expect_true(all(selected$selected_audit_min_coverage >= 0 &
                                selected$selected_audit_min_coverage <= 1))
    testthat::expect_true(all(selected$selected_ess_fraction > 0 &
                                selected$selected_ess_fraction <= 1))
  }
})

testthat::test_that("S32 confines the stress perturbation to invalid queries", {
  grid <- seq(30, 60, length.out = 4L)
  bank <- data.frame(z = grid)
  validity <- rbind(
    c(TRUE, TRUE, TRUE, TRUE),
    c(FALSE, TRUE, TRUE, TRUE),
    c(TRUE, TRUE, TRUE, FALSE),
    c(FALSE, TRUE, TRUE, TRUE),
    c(TRUE, TRUE, TRUE, FALSE)
  )
  predictions <- outer(seq(0.2, 0.4, length.out = 5L),
                       seq(0, 0.12, length.out = 4L), "+")
  attr(bank, "prediction_matrix") <- predictions
  attr(bank, "evaluation_validity") <- validity
  attr(bank, "observed_predictions") <- seq(0.2, 0.4, length.out = 5L)
  attr(bank, "observed_validity") <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
  attr(bank, "evaluation_response") <- c(0, 0, 1, 0, 1)
  attr(bank, "sample_focal") <- c(32, 38, 45, 52, 59)
  simulation <- run_s32(bank, logit_shift = 3, epsilon = 0.1)
  distances <- stats::setNames(
    simulation$sensitivity$centered_sup_distance,
    simulation$sensitivity$method
  )
  output <- tempfile(fileext = ".png")

  testthat::expect_gt(distances[["Ordinary PD"]], 0)
  testthat::expect_equal(distances[["Hard CSPD"]], 0, tolerance = 1e-12)
  testthat::expect_equal(distances[["PTPD"]], 0, tolerance = 1e-12)
  testthat::expect_gte(simulation$observed_unchanged_fraction, 0.6)
  testthat::expect_true(all(file.exists(plot_s32(simulation, output))))
})

testthat::test_that("S33 covariance path equals the projected-minus-original curve", {
  set.seed(33)
  predictions <- matrix(stats::rnorm(60), nrow = 12L, ncol = 5L)
  weights <- stats::rexp(12L)
  diagnostic <- s33_environment$empirical_reference_covariance(
    predictions, weights
  )

  testthat::expect_equal(
    diagnostic$covariance_path,
    diagnostic$direct_difference,
    tolerance = 1e-12
  )
  testthat::expect_lt(max(diagnostic$identity_error), 1e-12)
})

testthat::test_that("S34 reports distinct conditional-subgroup PDPs", {
  simulation <- run_s34(
    n_evaluation = 300L, grid = seq(-0.6, 0.6, length.out = 9L),
    seed = 34, minbucket_fraction = 0.10
  )
  output <- file.path(tempdir(), "s34-test")
  paths <- plot_s34(simulation, output)

  testthat::expect_gte(nrow(simulation$groups), 2L)
  testthat::expect_equal(sum(simulation$groups$n), 300L)
  testthat::expect_equal(sum(simulation$groups$fraction), 1, tolerance = 1e-12)
  testthat::expect_true(all(simulation$groups$x_q25 < simulation$groups$x_q75))
  testthat::expect_true(any(simulation$curves$dense))
  testthat::expect_true(all(file.exists(paths)))
})
