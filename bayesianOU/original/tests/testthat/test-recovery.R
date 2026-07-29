# Parameter-recovery regression test (Stan-based, SLOW).
#
# Disabled by default: it compiles the Stan model and runs MCMC. Enable with
#   Sys.setenv(BAYESOU_RUN_RECOVERY = "1")
# A heavier, harder version lives in validation/recovery_check.R.

test_that("the Stan model recovers known parameters (90% CI coverage)", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("BAYESOU_RUN_RECOVERY")),
              "Set BAYESOU_RUN_RECOVERY=1 to run the (slow) recovery test.")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if(check_stan_backend() == "none", "No Stan backend available.")

  set.seed(101)
  S <- 2L; T_obs <- 90L; T_train <- floor(T_obs * 0.70)

  truth <- list(
    theta_s = c(0.05, -0.05), kappa_s = c(0.30, 0.45),
    a3_s = -c(0.05, 0.07), beta0_s = c(0.25, 0.30),
    alpha_s = c(-0.5, -0.4), rho_s = c(0.80, 0.75),
    sigma_eta_s = c(0.30, 0.30), beta1 = 0.40, gamma = 0.30, nu = 6
  )

  ar1 <- function(n, phi) {
    x <- numeric(n); x[1] <- stats::rnorm(1, 0, 1 / sqrt(1 - phi^2))
    for (t in 2:n) x[t] <- phi * x[t - 1] + stats::rnorm(1)
    x
  }
  Xz <- sapply(seq_len(S), function(s) { x <- ar1(T_obs, 0.5); (x - mean(x)) / sd(x) })
  zTMG <- { x <- ar1(T_obs, 0.6); (x - mean(x)) / sd(x) }
  COM <- sapply(seq_len(S), function(s) pmax(1 + 0.4 * ar1(T_obs, 0.7), 0.05))
  K <- matrix(1, T_obs, S)
  cwm <- colMeans(COM[seq_len(T_train), , drop = FALSE])
  cws <- sapply(seq_len(S), function(s) sqrt(mean((COM[seq_len(T_train), s] - cwm[s])^2)))
  com_std <- sapply(seq_len(S), function(s) (COM[, s] - cwm[s]) / cws[s])

  h <- matrix(0, T_obs, S)
  for (s in seq_len(S)) {
    rho <- truth$rho_s[s]; hstd <- numeric(T_obs)
    hstd[1] <- stats::rnorm(1) / sqrt(1 - rho^2)
    for (t in 2:T_obs) hstd[t] <- rho * hstd[t - 1] + stats::rnorm(1)
    h[, s] <- truth$alpha_s[s] + truth$sigma_eta_s[s] * hstd
  }
  sig <- exp(0.5 * h)

  Yz <- matrix(0, T_obs, S); Yz[1, ] <- stats::rnorm(S)
  for (t in 2:T_obs) for (s in seq_len(S)) {
    zlag <- Yz[t - 1, s] - truth$theta_s[s]
    drift <- truth$kappa_s[s] * (truth$theta_s[s] - Yz[t - 1, s] + truth$a3_s[s] * zlag^3)
    betaT <- truth$beta0_s[s] + truth$beta1 * zTMG[t]
    mean_ <- drift + betaT * Xz[t - 1, s] + truth$gamma * com_std[t - 1, s]
    Yz[t, s] <- Yz[t - 1, s] + mean_ + sig[t, s] * stats::rt(1, truth$nu)
  }

  stan_dat <- list(
    n_levels = 1L,
    T = T_obs, S = S, T_train = T_train, T_lik = T_obs,
    Yz = Yz, Xz = Xz, zTMG_byK = zTMG, zTMG_exo = zTMG,
    soft_wedge = 0L, sigma_delta_z = 0.01, COM_ts = COM, K_ts = K,
    com_in_mean = 1L, mu_xz = rep(0, S),
    beta1_prior_mean = 0, beta1_prior_sd = 0.5,
    nu_prior_shape = 2, nu_prior_rate = 0.1,
    rho_prior_mean = 0.7, rho_prior_sd = 0.2,
    # Level-2 fields are length 0 in single-level mode.
    Phi_anchor_z = matrix(0.0, 0L, 0L), Gprime = numeric(0),
    theta_sep = 0L, sigma_phi_meas_prior_sd = 0.5,
    sigma_phi_meas_fixed = 0L, sigma_phi_meas_value = 0.5,
    V_anchor_z = matrix(0.0, 0L, 0L), kappa_cap = 2,
    # Per-level richness switches (canonical) and reconstruction fields (off).
    l1_cubic = 1L, l1_sv = 1L, l1_studentt = 1L, l1_hier = 1L,
    l2_cubic = 0L, l2_sv = 0L, l2_studentt = 0L, l2_hier = 1L,
    k_recon = 0L, k_cost = matrix(0.0, 0L, 0L), K_hat = matrix(0.0, 0L, 0L),
    Gprime_raw = numeric(0), phi_recon_center = numeric(0),
    phi_recon_scale = numeric(0), sigma_K_recon = 0.10
  )

  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(ou_nested_stan_code()),
    cpp_options = list(stan_threads = TRUE)
  )
  fit <- mod$sample(
    data = stan_dat, chains = 2, parallel_chains = 2, threads_per_chain = 2,
    iter_warmup = 600, iter_sampling = 600, adapt_delta = 0.9,
    max_treedepth = 10, seed = 123, refresh = 0
  )

  # Check only the structurally well-identified parameters. nu (Student-t df)
  # and the SV scale are weakly identified when stochastic volatility is present
  # (the prior dominates), so they are deliberately excluded from this
  # regression assertion to avoid flakiness; validation/recovery_check.R reports
  # the full picture including those nuisance parameters.
  truth_vec <- c(setNames(truth$kappa_s, sprintf("kappa_s[%d]", seq_len(S))),
                 beta1 = truth$beta1, gamma = truth$gamma)
  dm <- posterior::as_draws_matrix(fit$draws(names(truth_vec)))
  qs <- t(apply(dm[, names(truth_vec), drop = FALSE], 2, stats::quantile,
                probs = c(0.05, 0.95), na.rm = TRUE))
  covered <- truth_vec >= qs[, 1] & truth_vec <= qs[, 2]
  rh <- posterior::summarise_draws(fit$draws(names(truth_vec)), "rhat")$rhat

  expect_lt(max(rh, na.rm = TRUE), 1.1)              # chains mixed
  expect_gte(mean(covered), 0.6)                     # structural params covered

  # The reshape helper must also produce a coherent LOO over (T-1)*S obs.
  lr <- .compute_loo(fit, T_lik = T_obs, S = S)
  expect_s3_class(lr, "loo")
  expect_equal(nrow(lr$pointwise), (T_obs - 1L) * S)
})
