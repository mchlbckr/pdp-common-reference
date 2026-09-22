# Simulation S20: sensitivity of oracle support diagnostics to alpha.

s20_one <- function(alpha, repetition, n_evaluation, rho,
                    grid_half_width, grid_points, seed) {
  set.seed(seed + 100000L * as.integer(round(100 * alpha)) + repetition)
  evaluation <- simulate_s1_data(n_evaluation, rho)
  grid <- seq(-grid_half_width, grid_half_width, length.out = grid_points)
  radius <- sqrt(1 - rho^2) * stats::qnorm(1 - alpha / 2)
  oracle <- function(newdata) abs(newdata$x1 - rho * newdata$x2) <= radius
  h <- validity_matrix(evaluation, "x1", grid, oracle)
  isc <- colMeans(h)
  data.frame(
    alpha = alpha,
    repetition = repetition,
    mean_isc = mean(isc),
    min_isc = min(isc),
    max_isc = max(isc),
    gamma = mean(rowSums(h) == ncol(h))
  )
}

run_s20 <- function(alphas = c(0.01, 0.05, 0.10, 0.20), repetitions = 200L,
                    n_evaluation = 2000L, rho = 0.8,
                    grid_half_width = 1.1, grid_points = 21L,
                    seed = 20260925L,
                    workers = if (.Platform$OS.type == "windows") 1L else 8L) {
  design <- expand.grid(alpha = alphas, repetition = seq_len(repetitions))
  evaluate <- function(i) s20_one(
    design$alpha[[i]], design$repetition[[i]], n_evaluation, rho,
    grid_half_width, grid_points, seed
  )
  pieces <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(design)), evaluate, mc.cores = workers)
  } else lapply(seq_len(nrow(design)), evaluate)
  do.call(rbind, pieces)
}

summarize_s20 <- function(simulation) {
  measures <- c("mean_isc", "min_isc", "max_isc", "gamma")
  pieces <- split(simulation, simulation$alpha)
  do.call(rbind, lapply(pieces, function(part) {
    output <- data.frame(alpha = part$alpha[[1]])
    for (measure in measures) {
      output[[measure]] <- mean(part[[measure]])
      output[[paste0(measure, "_mc_se")]] <-
        stats::sd(part[[measure]]) / sqrt(nrow(part))
    }
    output
  }))
}
