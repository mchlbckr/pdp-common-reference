# Simulation S31: finite-sample certification of learned RCPD tilts.

run_s31 <- function(
  repetitions = 100L,
  audit_sizes = c(2500L, 5000L, 10000L, 20000L),
  bounds = c("empirical_bernstein", "hoeffding"),
  n_projection = 1000L, rho = 0.8, alpha = 0.05,
  fitted_epsilons = c(0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10),
  target_epsilon = 0.10, delta = 0.10,
  grid = seq(-1.1, 1.1, length.out = 21L),
  truth_size = 100000L, seed = 20261001L,
  workers = if (.Platform$OS.type == "windows") 1L else 8L
) {
  valid <- s14_oracle_validity(rho, alpha)
  truth_contexts <- data.frame(
    x1 = 0,
    x2 = stats::qnorm((seq_len(truth_size) - 0.5) / truth_size)
  )
  truth_validity <- validity_matrix(truth_contexts, "x1", grid, valid)
  design <- expand.grid(
    repetition = seq_len(repetitions),
    audit_size = audit_sizes,
    bound = bounds,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  maximum_audit_size <- max(audit_sizes)
  evaluate <- function(index) {
    row <- design[index, , drop = FALSE]
    # Generate one nested audit sample per repetition. Passing the maximum here
    # preserves pairing across all requested audit sizes.
    valid <- s14_oracle_validity(rho, alpha)
    projection <- simulate_s1_data(
      n_projection, rho, seed = seed + 100000L * row$repetition
    )
    audit_full <- simulate_s1_data(
      maximum_audit_size, rho,
      seed = seed + 200000L * row$repetition
    )
    audit <- audit_full[seq_len(row$audit_size), , drop = FALSE]
    projection_validity <- validity_matrix(projection, "x1", grid, valid)
    audit_validity <- validity_matrix(audit, "x1", grid, valid)
    certificate <- certify_calibrated_rcpd(
      projection_validity, audit_validity,
      fitted_epsilons = fitted_epsilons,
      target_epsilon = target_epsilon, delta = delta, bound = row$bound
    )
    naive <- fit_relaxed_reference(
      projection_validity, epsilon = target_epsilon, tolerance = 2e-5
    )
    naive_truth <- evaluate_rcpd_tilt(truth_validity, naive$multipliers)
    if (is.null(certificate$selected)) {
      return(data.frame(
        repetition = row$repetition, audit_size = row$audit_size,
        bound = row$bound, certified = FALSE,
        selected_epsilon = NA_real_, selected_projection_kl = NA_real_,
        selected_ess_fraction = NA_real_,
        selected_audit_min_coverage = NA_real_,
        selected_truth_min_coverage = NA_real_, certificate_holds = NA,
        naive_truth_min_coverage = naive_truth$min_coverage,
        naive_violation = naive_truth$min_coverage < 1 - target_epsilon,
        stringsAsFactors = FALSE
      ))
    }
    selected_truth <- evaluate_rcpd_tilt(
      truth_validity, certificate$selected$multipliers
    )
    data.frame(
      repetition = row$repetition, audit_size = row$audit_size,
      bound = row$bound, certified = TRUE,
      selected_epsilon = certificate$selected$fitted_epsilon,
      selected_projection_kl = certificate$selected$fit$kl,
      selected_ess_fraction =
        certificate$selected$fit$effective_sample_size / n_projection,
      selected_audit_min_coverage =
        certificate$selected$audit$audit_min_coverage,
      selected_truth_min_coverage = selected_truth$min_coverage,
      certificate_holds =
        selected_truth$min_coverage >= 1 - target_epsilon,
      naive_truth_min_coverage = naive_truth$min_coverage,
      naive_violation = naive_truth$min_coverage < 1 - target_epsilon,
      stringsAsFactors = FALSE
    )
  }
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(design)), evaluate, mc.cores = workers)
  } else {
    lapply(seq_len(nrow(design)), evaluate)
  }
  output <- do.call(rbind, pieces)
  attr(output, "settings") <- list(
    repetitions = repetitions, audit_sizes = audit_sizes, bounds = bounds,
    n_projection = n_projection, rho = rho, alpha = alpha,
    fitted_epsilons = fitted_epsilons,
    target_epsilon = target_epsilon, delta = delta,
    grid = grid, truth_size = truth_size, seed = seed
  )
  output
}

summarize_s31 <- function(simulation) {
  groups <- interaction(simulation$bound, simulation$audit_size, drop = TRUE)
  pieces <- split(simulation, groups)
  output <- do.call(rbind, lapply(pieces, function(part) {
    certified <- part[part$certified, , drop = FALSE]
    mean_or_na <- function(variable) {
      if (nrow(certified)) mean(certified[[variable]]) else NA_real_
    }
    mcse_or_na <- function(variable) {
      if (nrow(certified) > 1L) {
        stats::sd(certified[[variable]]) / sqrt(nrow(certified))
      } else {
        NA_real_
      }
    }
    data.frame(
      bound = part$bound[[1]], audit_size = part$audit_size[[1]],
      repetitions = nrow(part),
      certification_rate = mean(part$certified),
      certification_rate_mc_se =
        sqrt(mean(part$certified) * (1 - mean(part$certified)) / nrow(part)),
      violation_rate_among_certified = if (nrow(certified))
        mean(!certified$certificate_holds) else NA_real_,
      mean_selected_epsilon = mean_or_na("selected_epsilon"),
      mean_selected_epsilon_mc_se = mcse_or_na("selected_epsilon"),
      mean_selected_truth_coverage =
        mean_or_na("selected_truth_min_coverage"),
      mean_selected_truth_coverage_mc_se =
        mcse_or_na("selected_truth_min_coverage"),
      mean_selected_projection_kl = mean_or_na("selected_projection_kl"),
      mean_selected_ess_fraction = mean_or_na("selected_ess_fraction"),
      naive_violation_rate = mean(part$naive_violation),
      mean_naive_truth_coverage = mean(part$naive_truth_min_coverage),
      stringsAsFactors = FALSE
    )
  }))
  rownames(output) <- NULL
  output[order(output$bound, output$audit_size), , drop = FALSE]
}
