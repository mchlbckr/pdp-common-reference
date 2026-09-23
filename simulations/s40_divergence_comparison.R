# Simulation S40: which divergence should define the relaxed common reference?
#
# The relaxed reference is defined as the KL projection of the observed context
# distribution onto the query-validity constraint set. That choice is not
# forced. This study holds the constraint system fixed and varies only the
# divergence, comparing the KL projection with the chi-squared projection.
#
# One relation is analytic rather than empirical and must be stated before the
# results are read. With uniform base weights, the Pearson divergence from the
# base distribution is n * sum(w^2) - 1, and the effective sample size is
# 1 / sum(w^2). Minimising the former is therefore exactly maximising the
# latter, so the chi-squared projection attains the largest effective sample
# size of any feasible reference by construction. The open questions are how
# large the gap is at realistic operating points, and what the larger effective
# sample size costs in curve fidelity and in robustness to unsupported
# predictions.

s40_truth <- function(z) sin(pi * z / 2)

simulate_s40_data <- function(n, sigma = 0.45, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- stats::rnorm(n)
  data.frame(x = u + sigma * stats::rnorm(n), u = u)
}

s40_validity <- function(sigma = 0.45, alpha = 0.05) {
  radius <- stats::qnorm(1 - alpha / 2) * sigma
  function(newdata) abs(newdata$x - newdata$u) <= radius
}

s40_prediction <- function(direction = 0, sigma = 0.45, alpha = 0.05,
                           beta = 1.5, strength = 4) {
  radius <- stats::qnorm(1 - alpha / 2) * sigma
  function(newdata) {
    excess <- pmax(abs(newdata$x - newdata$u) - radius, 0)
    s40_truth(newdata$x) + beta * newdata$u +
      direction * strength * excess^2
  }
}

s40_center <- function(x) x - mean(x)

#' Fit both projections on one validity matrix and describe the weights.
s40_references <- function(h, epsilon, tolerance = 2e-7) {
  list(
    KL = fit_relaxed_reference(h, epsilon = epsilon, tolerance = tolerance),
    `Chi-squared` = fit_quadratic_reference(h, epsilon = epsilon,
                                            tolerance = tolerance)
  )
}

s40_one <- function(repetition, n_evaluation, sigma, alpha, epsilon, grid,
                    beta, strength, seed) {
  data <- simulate_s40_data(n_evaluation, sigma, seed + repetition)
  valid <- s40_validity(sigma, alpha)
  h <- validity_matrix(data, "x", grid, valid)
  epsilon_min <- minimum_relaxation_epsilon(h)
  if (epsilon + 1e-9 < epsilon_min) {
    return(NULL)
  }
  references <- s40_references(h, epsilon)
  truth_centered <- s40_center(s40_truth(grid))
  truth_amplitude <- diff(range(truth_centered))

  directions <- c(-1, 1)
  curves <- lapply(directions, function(direction) {
    prediction <- s40_prediction(direction, sigma, alpha, beta, strength)
    vapply(references, function(reference) {
      s40_center(estimate_weighted_pdp(
        data, "x", grid, prediction, reference$weights
      )$weighted_pdp)
    }, numeric(length(grid)))
  })

  do.call(rbind, lapply(names(references), function(name) {
    reference <- references[[name]]
    weights <- reference$weights
    centered <- curves[[2L]][, name]
    data.frame(
      repetition = repetition,
      epsilon = epsilon,
      sigma = sigma,
      divergence = name,
      ess_fraction = reference$effective_sample_size / n_evaluation,
      max_weight = max(weights) * n_evaluation,
      support_fraction = mean(weights > 0),
      kl_cost = reference$kl,
      chi_squared_cost = n_evaluation * sum(weights^2) - 1,
      min_coverage = reference$min_coverage,
      max_violation = reference$max_violation,
      normalized_shape_error = max(abs(centered - truth_centered)) /
        truth_amplitude,
      extension_sensitivity = max(abs(curves[[1L]][, name] -
                                        curves[[2L]][, name])),
      gamma = mean(rowSums(h) == ncol(h)),
      epsilon_min = epsilon_min,
      stringsAsFactors = FALSE
    )
  }))
}

run_s40 <- function(repetitions = 200L, n_evaluation = 1500L,
                    sigma_values = c(0.45, 0.8), alpha = 0.05,
                    epsilon_values = c(0.05, 0.10, 0.20, 0.30),
                    grid = seq(-0.75, 0.75, length.out = 21L),
                    beta = 1.5, strength = 4, seed = 20261203L) {
  designs <- expand.grid(
    sigma = sigma_values, epsilon = epsilon_values,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  results <- lapply(seq_len(nrow(designs)), function(index) {
    design <- designs[index, ]
    message(sprintf(
      "S40 design %d/%d (sigma %.2f, epsilon %.2f)",
      index, nrow(designs), design$sigma, design$epsilon
    ))
    parts <- lapply(seq_len(repetitions), function(repetition) {
      s40_one(
        repetition, n_evaluation, design$sigma, alpha, design$epsilon,
        grid, beta, strength, seed + 1000L * index
      )
    })
    do.call(rbind, parts)
  })
  do.call(rbind, results)
}

summarize_s40 <- function(simulation) {
  keys <- c("sigma", "epsilon", "divergence")
  metrics <- c("ess_fraction", "max_weight", "support_fraction", "kl_cost",
               "chi_squared_cost", "min_coverage", "normalized_shape_error",
               "extension_sensitivity", "gamma")
  grouped <- split(simulation, simulation[keys], drop = TRUE)
  summary <- do.call(rbind, lapply(grouped, function(part) {
    out <- part[1L, keys, drop = FALSE]
    out$repetitions <- nrow(part)
    for (metric in metrics) {
      out[[paste0("median_", metric)]] <- stats::median(part[[metric]])
    }
    out$mean_ess_fraction <- mean(part$ess_fraction)
    out$mean_ess_fraction_mc_se <-
      stats::sd(part$ess_fraction) / sqrt(nrow(part))
    out$mean_normalized_shape_error <- mean(part$normalized_shape_error)
    out$mean_normalized_shape_error_mc_se <-
      stats::sd(part$normalized_shape_error) / sqrt(nrow(part))
    out$mean_extension_sensitivity <- mean(part$extension_sensitivity)
    out$mean_extension_sensitivity_mc_se <-
      stats::sd(part$extension_sensitivity) / sqrt(nrow(part))
    out$max_constraint_violation <- max(part$max_violation)
    out
  }))
  rownames(summary) <- NULL
  summary <- summary[order(summary$sigma, summary$epsilon,
                           summary$divergence), ]

  # paired ESS gain of the chi-squared projection over the KL projection
  paired <- merge(
    simulation[simulation$divergence == "KL",
               c("repetition", "sigma", "epsilon", "ess_fraction",
                 "normalized_shape_error", "extension_sensitivity")],
    simulation[simulation$divergence == "Chi-squared",
               c("repetition", "sigma", "epsilon", "ess_fraction",
                 "normalized_shape_error", "extension_sensitivity")],
    by = c("repetition", "sigma", "epsilon"), suffixes = c("_kl", "_chi")
  )
  paired$ess_ratio <- paired$ess_fraction_chi / paired$ess_fraction_kl
  paired$shape_error_difference <-
    paired$normalized_shape_error_chi - paired$normalized_shape_error_kl
  paired$sensitivity_difference <-
    paired$extension_sensitivity_chi - paired$extension_sensitivity_kl
  contrast_groups <- split(paired, paired[c("sigma", "epsilon")], drop = TRUE)
  contrasts <- do.call(rbind, lapply(contrast_groups, function(part) {
    data.frame(
      sigma = part$sigma[[1L]],
      epsilon = part$epsilon[[1L]],
      repetitions = nrow(part),
      median_ess_ratio = stats::median(part$ess_ratio),
      mean_ess_ratio = mean(part$ess_ratio),
      mean_ess_ratio_mc_se = stats::sd(part$ess_ratio) / sqrt(nrow(part)),
      mean_shape_error_difference = mean(part$shape_error_difference),
      mean_shape_error_difference_mc_se =
        stats::sd(part$shape_error_difference) / sqrt(nrow(part)),
      mean_sensitivity_difference = mean(part$sensitivity_difference),
      mean_sensitivity_difference_mc_se =
        stats::sd(part$sensitivity_difference) / sqrt(nrow(part)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(contrasts) <- NULL

  list(summary = summary,
       contrasts = contrasts[order(contrasts$sigma, contrasts$epsilon), ])
}

#' Repeat the divergence comparison on the California Housing validity
#' matrices, which is the operating point where the weight concentration
#' reported in the application is most severe.
#'
#' The exploratory screening scripts own the data loader, the validity rules,
#' and the split rule, so they are loaded into a private environment exactly as
#' the confirmatory application does.
run_s40_california <- function(repetitions = 10L, base_seed = 20261101L,
                               epsilon_values = c(0.10, 0.15, 0.20),
                               alpha = 0.10, k = 21L) {
  environment_s40 <- new.env(parent = globalenv())
  for (script in c("s33_real_data_screen.R",
                   "s34_california_longitude_stability.R",
                   "s35_california_longitude_robustness.R")) {
    sys.source(file.path("exploration", script), envir = environment_s40)
  }
  data <- environment_s40$pilot_california_data()
  seeds <- base_seed + seq_len(repetitions) - 1L

  rows <- lapply(seq_along(seeds), function(index) {
    seed <- seeds[[index]]
    message(sprintf("S40 California split %d/%d", index, repetitions))
    split <- environment_s40$pilot_split(nrow(data$predictors), seed)
    split$evaluation <- environment_s40$pilot_cap_evaluation(
      split$evaluation, maximum = 5000L, seed = seed + 1L
    )
    training <- data$predictors[split$training, , drop = FALSE]
    calibration <- data$predictors[split$calibration, , drop = FALSE]
    evaluation <- data$predictors[split$evaluation, , drop = FALSE]
    limits <- stats::quantile(
      training$Longitude, c(0.40, 0.60), names = FALSE, type = 8
    )
    grid <- seq(limits[[1]], limits[[2]], length.out = k)
    learners <- list(
      "Locally scaled residual" = environment_s40$pilot_validity_rule(
        training, calibration, "Longitude", alpha
      ),
      "Conditional Gaussian density" =
        environment_s40$pilot_conditional_density_rule(
          training, calibration, "Longitude", alpha
        )
    )
    parts <- lapply(names(learners), function(learner) {
      validity <- validity_matrix(
        evaluation, "Longitude", grid, learners[[learner]]
      )
      n <- nrow(validity)
      epsilon_min <- minimum_relaxation_epsilon(validity)
      usable <- epsilon_values[epsilon_values + 1e-9 >= epsilon_min]
      if (length(usable) == 0L) {
        return(NULL)
      }
      do.call(rbind, lapply(usable, function(epsilon) {
        references <- list(
          KL = environment_s40$fit_relaxed_reference_precise(
            validity, epsilon = epsilon
          ),
          `Chi-squared` = fit_quadratic_reference(validity, epsilon = epsilon)
        )
        do.call(rbind, lapply(names(references), function(name) {
          reference <- references[[name]]
          weights <- reference$weights
          data.frame(
            split_id = index, seed = seed, learner = learner,
            epsilon = epsilon, divergence = name,
            n_evaluation = n,
            ess_fraction = reference$effective_sample_size / n,
            max_weight = max(weights) * n,
            support_fraction = mean(weights > 0),
            kl_cost = reference$kl,
            chi_squared_cost = n * sum(weights^2) - 1,
            min_coverage = min(as.numeric(crossprod(weights,
                                                    unclass(validity) * 1))),
            simultaneous_valid_mass = mean(rowSums(validity) == ncol(validity)),
            stringsAsFactors = FALSE
          )
        }))
      }))
    })
    do.call(rbind, parts)
  })
  do.call(rbind, rows)
}

summarize_s40_california <- function(simulation) {
  keys <- c("learner", "epsilon", "divergence")
  metrics <- c("ess_fraction", "max_weight", "support_fraction", "kl_cost",
               "chi_squared_cost", "min_coverage", "simultaneous_valid_mass")
  grouped <- split(simulation, simulation[keys], drop = TRUE)
  summary <- do.call(rbind, lapply(grouped, function(part) {
    out <- part[1L, keys, drop = FALSE]
    out$splits <- nrow(part)
    for (metric in metrics) {
      out[[paste0("median_", metric)]] <- stats::median(part[[metric]])
    }
    out
  }))
  rownames(summary) <- NULL

  paired <- merge(
    simulation[simulation$divergence == "KL",
               c("split_id", "learner", "epsilon", "ess_fraction")],
    simulation[simulation$divergence == "Chi-squared",
               c("split_id", "learner", "epsilon", "ess_fraction")],
    by = c("split_id", "learner", "epsilon"), suffixes = c("_kl", "_chi")
  )
  paired$ess_ratio <- paired$ess_fraction_chi / paired$ess_fraction_kl
  contrast_groups <- split(paired, paired[c("learner", "epsilon")], drop = TRUE)
  contrasts <- do.call(rbind, lapply(contrast_groups, function(part) {
    data.frame(
      learner = part$learner[[1L]], epsilon = part$epsilon[[1L]],
      splits = nrow(part),
      median_ess_ratio = stats::median(part$ess_ratio),
      mean_ess_ratio = mean(part$ess_ratio),
      mean_ess_ratio_mc_se = stats::sd(part$ess_ratio) / sqrt(nrow(part)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(contrasts) <- NULL
  list(summary = summary[order(summary$learner, summary$epsilon,
                               summary$divergence), ],
       contrasts = contrasts[order(contrasts$learner, contrasts$epsilon), ])
}
