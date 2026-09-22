source(file.path("..", "..", "R", "support.R"))
source(file.path("..", "..", "R", "validity.R"))

testthat::test_that("split conformal validity returns a logical decision per row", {
  set.seed(1)
  data <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
  model <- fit_split_conformal_validity(data[1:60, ], data[61:90, ], "x1")
  accepted <- predict_validity(model, data[91:120, ])

  testthat::expect_type(accepted, "logical")
  testthat::expect_length(accepted, 30)
  testthat::expect_false(anyNA(accepted))
  testthat::expect_lte(model$calibration_exceedance, model$alpha)
})

testthat::test_that("conditional density validity uses a calibrated score threshold", {
  set.seed(2)
  data <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
  model <- fit_conditional_density_validity(data[1:60, ], data[61:90, ], "x1")
  accepted <- predict_conditional_density_validity(model, data[91:120, ])

  testthat::expect_s3_class(model, "conditional_density_validity")
  testthat::expect_type(accepted, "logical")
  testthat::expect_length(accepted, 30)
  testthat::expect_false(anyNA(accepted))
  testthat::expect_lte(model$calibration_exceedance, model$alpha)
})

testthat::test_that("joint Gaussian validity detects globally unusual contexts", {
  set.seed(101)
  training <- data.frame(x = rnorm(200), u = rnorm(200))
  calibration <- data.frame(x = rnorm(100), u = rnorm(100))
  fit <- fit_joint_gaussian_validity(training, calibration, alpha = 0.1)
  ordinary <- data.frame(x = 0, u = 0)
  unusual <- data.frame(x = 8, u = 8)

  testthat::expect_identical(predict_joint_gaussian_validity(fit, ordinary), TRUE)
  testthat::expect_identical(predict_joint_gaussian_validity(fit, unusual), FALSE)
  testthat::expect_lte(fit$calibration_exceedance, fit$alpha)
})

testthat::test_that("locally scaled residual validity adapts its width", {
  set.seed(12)
  training <- data.frame(x = rnorm(120), u = runif(120, -1, 1))
  training$x <- training$u + exp(training$u) * training$x
  calibration <- data.frame(x = rnorm(100), u = runif(100, -1, 1))
  calibration$x <- calibration$u + exp(calibration$u) * calibration$x
  fit <- fit_locally_scaled_validity(training, calibration, "x", alpha = 0.1)
  probes <- data.frame(x = c(0, 0), u = c(-1, 1))

  testthat::expect_s3_class(fit, "locally_scaled_validity")
  testthat::expect_length(predict_locally_scaled_validity(fit, probes), 2)
  testthat::expect_true(is.finite(fit$cutoff))
  testthat::expect_lte(fit$calibration_exceedance, fit$alpha)
})

testthat::test_that("kNN validity can calibrate on an independent split", {
  training <- data.frame(x = seq(0, 1, length.out = 20), u = seq(0, 1, length.out = 20))
  calibration <- data.frame(x = seq(0.025, 0.975, length.out = 18), u = seq(0.025, 0.975, length.out = 18))
  fit <- fit_knn_validity(training, k = 3, alpha = 0.1, calibration_data = calibration)

  testthat::expect_s3_class(fit, "knn_validity")
  testthat::expect_true(is.finite(fit$threshold))
  testthat::expect_lte(fit$calibration_exceedance, 0.1)
})
