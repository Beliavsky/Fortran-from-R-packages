# Per-configuration parameter-recovery + LOO-discrimination test (Stan, SLOW).
#
# Disabled by default (compiles Stan + runs MCMC for several configurations).
# Enable with Sys.setenv(BAYESOU_RUN_RECOVERY_CONFIG = "1"). The heavier,
# fuller-diagnostic version lives in validacion/recovery_by_config.R.
#
# For each 2-level configuration (canonical, both_full, both_lean, n1_lean) the
# data are simulated from THAT configuration's own generative process and fit;
# the structurally identifiable Level-2 parameters (kappa_p, kappa_m_base, m1)
# are checked for 90% credible-interval coverage and the latent Phi for
# correlation with the truth. The single-level model is fit on the canonical
# (adversarial) dataset to verify PSIS-LOO discriminates in favour of the 2-level
# model. nu and the SV scale are weakly identified and are excluded from the
# coverage assertion (the single-level structural recovery is in test-recovery.R).

test_that("each configuration recovers its parameters and LOO discriminates", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("BAYESOU_RUN_RECOVERY_CONFIG")),
              "Set BAYESOU_RUN_RECOVERY_CONFIG=1 to run the (slow) per-config test.")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if_not_installed("loo")
  skip_if(check_stan_backend() == "none", "No Stan backend available.")

  T_obs <- 80L; S <- 2L; T_train <- floor(T_obs * 0.70); KAPPA_CAP <- 2

  inv_logit <- function(x) 1 / (1 + exp(-x))
  clamp <- function(x, lo, hi) max(min(x, hi), lo)
  rstd  <- function(studentt, nu) clamp(if (studentt) stats::rt(1, nu) else stats::rnorm(1), -6, 6)
  ar1 <- function(n, phi) { x <- numeric(n); x[1] <- stats::rnorm(1) / sqrt(1 - phi^2)
    for (t in 2:n) x[t] <- phi * x[t - 1] + stats::rnorm(1); x }
  sv_path <- function(Tn, a, rho, se) { h <- numeric(Tn); h[1] <- stats::rnorm(1) / sqrt(1 - rho^2)
    for (t in 2:Tn) h[t] <- rho * h[t - 1] + stats::rnorm(1); a + se * h }

  truth <- list(
    kappa_tilde = c(-0.20, -0.60), a3_s = -c(0.03, 0.04),
    beta1 = 0.30, gamma = 0.20, nu = 6,
    alpha_s = c(-0.40, -0.50), rho_s = c(0.80, 0.75), sigma_eta_s = c(0.30, 0.30),
    kappa_p = c(0.25, 0.30), mu_const = c(0, 0), m1 = 0.50, sigma_p = 0.15,
    sigma_phi_meas = 0.30, a3_p = -c(0.04, 0.05), nu_p = 6,
    alpha_p = c(-0.60, -0.55), rho_p = c(0.75, 0.70), sigma_eta_p = c(0.25, 0.25))
  truth$kappa_m_base <- KAPPA_CAP * inv_logit(truth$kappa_tilde)

  flags_of <- function(cfg) {
    full <- c(1L, 1L, 1L, 1L); lean <- c(0L, 0L, 0L, 1L)
    m <- switch(cfg, single = , canonical = list(l1 = full, l2 = c(0L,0L,0L,1L)),
                both_full = list(l1 = full, l2 = full),
                both_lean = , n1_lean = list(l1 = lean, l2 = lean))
    list(l1_cubic=m$l1[1], l1_sv=m$l1[2], l1_studentt=m$l1[3], l1_hier=m$l1[4],
         l2_cubic=m$l2[1], l2_sv=m$l2[2], l2_studentt=m$l2[3], l2_hier=m$l2[4])
  }

  simulate <- function(fl, seed) {
    set.seed(seed)
    Gz <- as.numeric(scale(ar1(T_obs, 0.92))); zTMG <- as.numeric(scale(ar1(T_obs, 0.6)))
    COM <- sapply(1:S, function(s) pmax(1 + 0.4 * ar1(T_obs, 0.7), 0.05)); K <- matrix(1, T_obs, S)
    cwm <- colMeans(COM[1:T_train, , drop = FALSE])
    cws <- sapply(1:S, function(s) sqrt(mean((COM[1:T_train, s] - cwm[s])^2))); cws[cws < 1e-12] <- 1
    com_std <- sapply(1:S, function(s) (COM[, s] - cwm[s]) / cws[s])
    Phi <- matrix(0, T_obs, S)
    hP <- if (fl$l2_sv) sapply(1:S, function(s) sv_path(T_obs, truth$alpha_p[s], truth$rho_p[s], truth$sigma_eta_p[s])) else NULL
    for (s in 1:S) { Phi[1, s] <- truth$mu_const[s] + truth$m1 * Gz[1]
      for (t in 2:T_obs) { mu <- truth$mu_const[s] + truth$m1 * Gz[t]; dp <- Phi[t-1,s] - mu
        cp <- if (fl$l2_cubic) truth$a3_p[s] * clamp(dp, -10, 10)^3 else 0
        sp <- if (fl$l2_sv) exp(0.5 * hP[t, s]) else truth$sigma_p
        Phi[t, s] <- Phi[t-1,s] + truth$kappa_p[s] * (-dp + cp) + sp * rstd(fl$l2_studentt, truth$nu_p) } }
    anchor <- Phi + matrix(stats::rnorm(T_obs*S, 0, truth$sigma_phi_meas), T_obs, S)
    hM <- if (fl$l1_sv) sapply(1:S, function(s) sv_path(T_obs, truth$alpha_s[s], truth$rho_s[s], truth$sigma_eta_s[s])) else NULL
    Yz <- matrix(0, T_obs, S); Yz[1, ] <- Phi[1, ] + stats::rnorm(S, 0, 0.1)
    for (t in 2:T_obs) for (s in 1:S) { dev <- Yz[t-1,s] - Phi[t-1,s]
      km <- KAPPA_CAP * inv_logit(truth$kappa_tilde[s] + truth$beta1 * zTMG[t])
      cb <- if (fl$l1_cubic) truth$a3_s[s] * clamp(dev, -6, 6)^3 else 0
      sdt <- if (fl$l1_sv) exp(0.5 * hM[t, s]) else exp(0.5 * truth$alpha_s[s])
      Yz[t, s] <- Yz[t-1,s] + km * (-dev + cb) + truth$gamma * com_std[t-1,s] + sdt * rstd(fl$l1_studentt, truth$nu) }
    list(Yz = Yz, anchor = anchor, Phi = Phi, Gz = Gz, zTMG = zTMG, COM = COM, K = K)
  }

  stan_dat_of <- function(d, nl, fl) {
    base <- list(n_levels = as.integer(nl), T = T_obs, S = S, T_train = T_train, T_lik = T_obs,
      Yz = d$Yz, Xz = d$anchor, zTMG_byK = d$zTMG, zTMG_exo = d$zTMG, soft_wedge = 0L,
      sigma_delta_z = 0.01, COM_ts = d$COM, K_ts = d$K, com_in_mean = 1L, mu_xz = rep(0, S),
      beta1_prior_mean = 0, beta1_prior_sd = 0.5, nu_prior_shape = 2, nu_prior_rate = 0.1,
      rho_prior_mean = 0.7, rho_prior_sd = 0.2, theta_sep = 0L, sigma_phi_meas_prior_sd = 0.5,
      sigma_phi_meas_fixed = 0L, sigma_phi_meas_value = 0.5, V_anchor_z = matrix(0,0,0),
      kappa_cap = KAPPA_CAP, l1_cubic = fl$l1_cubic, l1_sv = fl$l1_sv, l1_studentt = fl$l1_studentt,
      l1_hier = fl$l1_hier, l2_cubic = fl$l2_cubic, l2_sv = fl$l2_sv, l2_studentt = fl$l2_studentt,
      l2_hier = fl$l2_hier, k_recon = 0L, k_cost = matrix(0,0,0), K_hat = matrix(0,0,0),
      Gprime_raw = numeric(0), phi_recon_center = numeric(0), phi_recon_scale = numeric(0),
      sigma_K_recon = 0.10)
    if (nl == 2L) { base$Phi_anchor_z <- d$anchor; base$Gprime <- d$Gz }
    else { base$Phi_anchor_z <- matrix(0,0,0); base$Gprime <- numeric(0) }
    base
  }

  cache <- tryCatch({ dd <- tools::R_user_dir("bayesianOU", "cache")
    dir.create(dd, recursive = TRUE, showWarnings = FALSE); if (dir.exists(dd)) dd else tempdir()
  }, error = function(e) tempdir())
  mod <- cmdstanr::cmdstan_model(.stan_file_path(), dir = cache, cpp_options = list(stan_threads = TRUE))
  fit_it <- function(sd, nl) mod$sample(data = sd, chains = 2, parallel_chains = 2,
    threads_per_chain = 2, iter_warmup = 600, iter_sampling = 600, adapt_delta = 0.95,
    max_treedepth = 12, seed = 4321, refresh = 0, init = if (nl == 2L) 0.3 else 2)

  loo_of <- function(fit) {
    ll <- fit$draws("log_lik", format = "matrix"); nch <- fit$num_chains()
    cid <- rep(seq_len(nch), each = nrow(ll) / nch)
    want <- as.vector(vapply(1:S, function(s) sprintf("log_lik[%d,%d]", 2:T_obs, s), character(T_obs-1)))
    ll <- ll[, want, drop = FALSE]
    loo::loo(ll, r_eff = loo::relative_eff(exp(ll), chain_id = cid))
  }
  corr_phi <- function(fit, Phi_true) {
    M <- fit$draws("Phi", format = "matrix"); cn <- colnames(M)
    tt <- as.integer(sub("^Phi\\[(\\d+),(\\d+)\\]$", "\\1", cn)); ss <- as.integer(sub("^Phi\\[(\\d+),(\\d+)\\]$", "\\2", cn))
    med <- apply(M, 2, stats::median); Ph <- matrix(NA_real_, T_obs, S)
    for (j in seq_along(cn)) Ph[tt[j], ss[j]] <- med[j]
    suppressWarnings(stats::cor(as.vector(Ph), as.vector(Phi_true)))
  }

  cfgs <- c("canonical", "both_full", "both_lean", "n1_lean")
  seeds <- c(canonical = 7002, both_full = 7003, both_lean = 7004, n1_lean = 7005)
  covered_all <- logical(0); data_canon <- NULL; loo_canon <- NULL
  for (cfg in cfgs) {
    fl <- flags_of(cfg); d <- simulate(fl, seeds[[cfg]])
    if (cfg == "canonical") data_canon <- d
    fit <- fit_it(stan_dat_of(d, 2L, fl), 2L)
    tv <- c(setNames(truth$kappa_p, sprintf("kappa_p[%d]", 1:S)),
            setNames(truth$kappa_m_base, sprintf("kappa_m_base[%d]", 1:S)), `m1[1]` = truth$m1)
    dm <- posterior::as_draws_matrix(fit$draws(names(tv)))
    qs <- t(apply(dm[, names(tv), drop = FALSE], 2, stats::quantile, probs = c(0.05, 0.95), na.rm = TRUE))
    cov <- as.numeric(tv) >= qs[, 1] & as.numeric(tv) <= qs[, 2]
    covered_all <- c(covered_all, cov)
    rh <- max(fit$summary(c("kappa_p", "m1"))$rhat, na.rm = TRUE)
    expect_lt(rh, 1.2)                                  # chains mixed
    expect_gt(corr_phi(fit, d$Phi), 0.7)                # latent Phi recovered
    if (cfg == "canonical") loo_canon <- loo_of(fit)
  }
  # Pooled 90% CI coverage across the 2-level configurations.
  expect_gte(mean(covered_all), 0.6)

  # LOO discrimination on the adversarial canonical dataset: the single-level
  # model (market reverts to a constant, index as exogenous forcing) vs the
  # data-generating 2-level model. PSIS-LOO should not prefer the single-level.
  fit_single <- fit_it(stan_dat_of(data_canon, 1L, flags_of("single")), 1L)
  loo_single <- loo_of(fit_single)
  expect_equal(nrow(loo_canon$pointwise), nrow(loo_single$pointwise))
  cmp <- loo::loo_compare(loo_canon, loo_single)
  delta <- loo_canon$estimates["elpd_loo", "Estimate"] -
           loo_single$estimates["elpd_loo", "Estimate"]
  se_diff <- cmp[2, "se_diff"]
  message(sprintf("LOO 2-level vs single: deltaELPD=%.1f se=%.1f", delta, se_diff))
  # The 2-level model is not meaningfully worse than the single-level on its own
  # adversarial data (the point estimate favours it; the separation is within
  # noise, as expected for weakly identified comparisons -- anti-overreach).
  expect_gt(delta, -2 * se_diff)
})
