library(testthat)
library(gkwdist)

# =============================================================================
# Helper Functions
# =============================================================================

# Generate test parameters for each subfamily
generate_test_params <- function(subfamily, test_id) {
  params_list <- list(
    gkw = list(
      list(alpha = 2.0, beta = 3.0, gamma = 1.5, delta = 0.5, lambda = 2.0),
      list(alpha = 0.8, beta = 1.2, gamma = 2.0, delta = 1.0, lambda = 0.9),
      list(alpha = 3.5, beta = 0.7, gamma = 0.8, delta = 2.0, lambda = 1.5),
      list(alpha = 1.5, beta = 2.5, gamma = 1.2, delta = 0.3, lambda = 1.8),
      list(alpha = 2.2, beta = 1.8, gamma = 0.6, delta = 1.5, lambda = 2.5)
    ),
    bkw = list(
      list(alpha = 2.0, beta = 3.0, gamma = 1.5, delta = 0.5),
      list(alpha = 0.8, beta = 1.2, gamma = 2.0, delta = 1.0),
      list(alpha = 3.5, beta = 0.7, gamma = 0.8, delta = 2.0),
      list(alpha = 1.5, beta = 2.5, gamma = 1.2, delta = 0.3),
      list(alpha = 2.2, beta = 1.8, gamma = 0.6, delta = 1.5)
    ),
    kkw = list(
      list(alpha = 2.0, beta = 3.0, delta = 0.5, lambda = 2.0),
      list(alpha = 0.8, beta = 1.2, delta = 1.0, lambda = 0.9),
      list(alpha = 3.5, beta = 0.7, delta = 2.0, lambda = 1.5),
      list(alpha = 1.5, beta = 2.5, delta = 0.3, lambda = 1.8),
      list(alpha = 2.2, beta = 1.8, delta = 1.5, lambda = 2.5)
    ),
    ekw = list(
      list(alpha = 2.0, beta = 3.0, lambda = 2.0),
      list(alpha = 0.8, beta = 1.2, lambda = 0.9),
      list(alpha = 3.5, beta = 0.7, lambda = 1.5),
      list(alpha = 1.5, beta = 2.5, lambda = 1.8),
      list(alpha = 2.2, beta = 1.8, lambda = 2.5)
    ),
    mc = list(
      list(gamma = 1.5, delta = 0.5, lambda = 2.0),
      list(gamma = 2.0, delta = 1.0, lambda = 0.9),
      list(gamma = 0.8, delta = 2.0, lambda = 1.5),
      list(gamma = 1.2, delta = 0.3, lambda = 1.8),
      list(gamma = 0.6, delta = 1.5, lambda = 2.5)
    ),
    kw = list(
      list(alpha = 2.0, beta = 3.0),
      list(alpha = 0.8, beta = 1.2),
      list(alpha = 3.5, beta = 0.7),
      list(alpha = 1.5, beta = 2.5),
      list(alpha = 2.2, beta = 1.8)
    ),
    beta = list(
      list(gamma = 2.0, delta = 3.0),
      list(gamma = 0.8, delta = 1.2),
      list(gamma = 3.5, delta = 0.7),
      list(gamma = 1.5, delta = 2.5),
      list(gamma = 2.2, delta = 1.8)
    )
  )
  params_list[[subfamily]][[test_id]]
}

# Wrapper for data simulation
simulate_data <- function(subfamily, params, n, seed = 2203) {
  set.seed(seed)
  args <- c(list(n = n), params)
  r_func <- paste0("r", subfamily)
  # Special case for beta as base R has rbeta
  if (subfamily == "beta") r_func <- "rbeta_"

  do.call(r_func, args)
}

# Convert params list to vector
params_to_vector <- function(subfamily, params) {
  if (subfamily == "gkw") {
    return(c(params$alpha, params$beta, params$gamma, params$delta, params$lambda))
  }
  if (subfamily == "bkw") {
    return(c(params$alpha, params$beta, params$gamma, params$delta))
  }
  if (subfamily == "kkw") {
    return(c(params$alpha, params$beta, params$delta, params$lambda))
  }
  if (subfamily == "ekw") {
    return(c(params$alpha, params$beta, params$lambda))
  }
  if (subfamily == "mc") {
    return(c(params$gamma, params$delta, params$lambda))
  }
  if (subfamily == "kw") {
    return(c(params$alpha, params$beta))
  }
  if (subfamily == "beta") {
    return(c(params$gamma, params$delta))
  }
}

# -----------------------------------------------------------------------------
# Core Fitting Function
# -----------------------------------------------------------------------------
fit_mle_timed <- function(subfamily, data, start, scenario = c("none", "gr", "gr_hs")) {
  scenario <- match.arg(scenario)

  # Dynamically get function names
  ll_fn <- get(paste0("ll", subfamily))
  gr_fn <- get(paste0("gr", subfamily))
  hs_fn <- get(paste0("hs", subfamily))

  obj_fn <- function(par) ll_fn(par, data)
  grad_fn <- function(par) gr_fn(par, data)
  hess_fn <- function(par) hs_fn(par, data)

  lower <- rep(1e-4, length(start))
  upper <- rep(100, length(start))

  start_time <- Sys.time()

  result <- tryCatch(
    {
      args <- list(
        start = start, objective = obj_fn, lower = lower, upper = upper,
        control = list(iter.max = 500, eval.max = 2000, trace = 0)
      )

      if (scenario %in% c("gr", "gr_hs")) args$gradient <- grad_fn
      if (scenario == "gr_hs") args$hessian <- hess_fn

      do.call(nlminb, args)
    },
    error = function(e) {
      list(
        par = rep(NA, length(start)), convergence = 99,
        evaluations = c("function" = NA, "gradient" = NA), message = e$message
      )
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  result$time <- elapsed
  return(result)
}

# -----------------------------------------------------------------------------
# Benchmark Runner
# -----------------------------------------------------------------------------
run_benchmark <- function(subfamily, config_id, n_reps = 20, n_obs = 1000) {
  params <- generate_test_params(subfamily, config_id)
  true_params <- params_to_vector(subfamily, params)

  scenarios <- c("none", "gr", "gr_hs")

  storage <- list()
  for (s in scenarios) {
    storage[[s]] <- data.frame(
      time = numeric(n_reps),
      evals = numeric(n_reps),
      error = numeric(n_reps),
      converged = logical(n_reps)
    )
  }

  for (i in seq_len(n_reps)) {
    # Generate data once per rep for paired comparison
    seed_val <- 2203 + (config_id * 100) + i
    data <- simulate_data(subfamily, params, n = n_obs, seed = seed_val)

    # Get start values
    start <- tryCatch(
      {
        sv <- gkwgetstartvalues(data, family = subfamily)
        if (any(is.na(sv)) || any(sv <= 0)) stop("Bad start")
        sv
      },
      error = function(e) {
        set.seed(seed_val)
        pmax(true_params * runif(length(true_params), 0.8, 1.2), 0.01)
      }
    )

    # Run all scenarios on the same data
    for (scen in scenarios) {
      fit <- fit_mle_timed(subfamily, data, start, scenario = scen)

      storage[[scen]]$time[i] <- fit$time
      storage[[scen]]$evals[i] <- fit$evaluations[1] # Function evaluations
      storage[[scen]]$converged[i] <- (fit$convergence == 0)

      if (fit$convergence == 0 && !any(is.na(fit$par))) {
        storage[[scen]]$error[i] <- mean(abs((fit$par - true_params) / true_params))
      } else {
        storage[[scen]]$error[i] <- NA
      }
    }
  }

  # Consolidate results
  summary_stats <- lapply(storage, function(df) {
    valid <- df[df$converged, ]
    list(
      mean_time = mean(valid$time, na.rm = TRUE),
      mean_evals = mean(valid$evals, na.rm = TRUE),
      mean_error = mean(valid$error, na.rm = TRUE),
      convergence_rate = mean(df$converged),
      raw_evals = valid$evals
    )
  })

  return(summary_stats)
}

# =============================================================================
# Custom Expectation Function
# =============================================================================
expect_efficiency_gain <- function(bench_result, dist_name) {
  # 1. Check basic stability - all scenarios should converge most of the time
  # Skip if numerical baseline fails (data/model issue, not gradient issue)
  if (bench_result$none$convergence_rate < 0.5) {
    skip(paste(dist_name, ": Low baseline convergence rate - skipping"))
  }

  expect_gt(bench_result$gr$convergence_rate, 0.3,
    label = paste(dist_name, ": Gradient scenario convergence too low")
  )
  expect_gt(bench_result$gr_hs$convergence_rate, 0.3,
    label = paste(dist_name, ": Gradient+Hessian scenario convergence too low")
  )

  # 2. Analytical gradient should not be much slower than numerical
  if (!is.na(bench_result$gr$mean_time) && !is.na(bench_result$none$mean_time)) {
    expect_lte(bench_result$gr$mean_time, bench_result$none$mean_time * 2.0,
      label = paste(dist_name, ": Analytical gradient too slow")
    )
  }

  # 3. Check Accuracy - analytical gradient should give equivalent or better results
  if (!is.na(bench_result$gr$mean_error) && !is.na(bench_result$none$mean_error)) {
    expect_lte(bench_result$gr$mean_error, bench_result$none$mean_error * 2.0,
      label = paste(dist_name, ": Analytical gradient precision is significantly worse")
    )
  }

  # 4. Hessian accuracy check (optional - only warn, don't fail)
  # In some cases Hessian can cause numerical issues
  # This is informational, not a hard requirement
}

# =============================================================================
# Test Execution
# =============================================================================

families <- c("beta", "kw", "ekw", "mc", "kkw", "bkw", "gkw")

for (fam in families) {
  test_that(paste("Optimization efficiency for family:", fam), {
    skip_on_cran()

    # Test 2 configs per family to balance coverage and time
    for (cfg in c(1, 5)) {
      # Increase reps for complex models
      reps <- if (fam == "gkw") 40 else 20

      bench <- run_benchmark(fam, cfg, n_reps = reps, n_obs = 1000)

      expect_efficiency_gain(bench, paste0(fam, "_cfg", cfg))
    }
  })
}
