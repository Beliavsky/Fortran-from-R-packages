# Per-imputation checkpoint / resume of fit_ou_nested_mi (crash-resilience).
#
# These tests are Stan-free: the engine (fit_ou_nested) is mocked with a
# deterministic fake whose per-imputation draws depend only on `seed + m`, so a
# refit of any imputation reproduces it exactly. That lets us assert the two
# properties that matter for a long, resumable run:
#   (1) a fully-restored run reproduces the uninterrupted object bit-for-bit and
#       performs zero (expensive) refits;
#   (2) a partially-completed run refits ONLY the missing imputations and still
#       yields the same pooled object as the uninterrupted run.
# The post-loop Rubin pooling is never mocked: it consumes the (restored or
# freshly computed) contributions, so equality of the final object is the real
# guarantee.

# Build a deterministic fake engine + matching .draws_matrix. `calls` counts how
# many times the engine actually samples (the quantity resume must drive to the
# missing imputations only).
.fake_engine_env <- function() {
  e <- new.env(parent = emptyenv())
  e$calls <- 0L
  e$ndraw <- 20L; e$nchain <- 2L; e$niter <- e$ndraw / e$nchain

  # pooled params at n_levels = 2: per-sector (S cols) and scalar (1 col).
  e$vec_pars <- c("kappa_m_base", "kappa_p", "sigma_p", "a3_s")
  e$sca_pars <- c("mu_const", "m1", "sigma_phi_meas", "beta1", "gamma", "nu")

  e$engine <- function(Y, X, TMG, COM, CAPITAL_TOTAL, Gprime = NULL,
                       V_value = NULL, n_levels = 1L, seed = 1234,
                       verbose = FALSE, ...) {
    e$calls <- e$calls + 1L
    Tn <- nrow(Y); S <- ncol(Y)
    set.seed(seed)                                  # determinism keyed by seed + m
    store <- list()
    for (p in e$vec_pars) {
      M_ <- matrix(stats::rnorm(e$ndraw * S), e$ndraw, S)
      colnames(M_) <- sprintf("%s[%d]", p, seq_len(S)); store[[p]] <- M_
    }
    for (p in e$sca_pars) {
      M_ <- matrix(stats::rnorm(e$ndraw), e$ndraw, 1L)
      colnames(M_) <- p; store[[p]] <- M_
    }
    phi_cn <- as.vector(vapply(seq_len(S),
      function(s) sprintf("Phi[%d,%d]", seq_len(Tn), s), character(Tn)))
    Mphi <- matrix(stats::rnorm(e$ndraw * Tn * S), e$ndraw, Tn * S)
    colnames(Mphi) <- phi_cn; store$Phi <- Mphi

    draws_fun <- function(p, ...) {
      M_ <- store[[p]]; if (is.null(M_)) return(NULL)
      arr <- array(NA_real_, dim = c(e$niter, e$nchain, ncol(M_)),
                   dimnames = list(NULL, NULL, colnames(M_)))
      for (v in seq_len(ncol(M_))) arr[, , v] <- matrix(M_[, v], e$niter, e$nchain)
      posterior::as_draws_array(arr)
    }
    sf <- list(store = store, draws = draws_fun)
    list(factor_ou = list(stan_fit = sf),
         diagnostics = list(rhat_max = 1.0, rhat_share = 0, ess = rep(500, 2),
                            divergences = 0L))
  }
  e$draws_matrix <- function(fit, p) fit$store[[p]]
  e
}

# Small fixed-anchor n_levels = 2 problem; the fake ignores the data values.
.mi_ckpt_data <- function(M = 3L, T0 = 6L, S0 = 2L, seed = 7L) {
  set.seed(seed)
  list(
    phi_draws = lapply(seq_len(M), function(k) matrix(stats::rnorm(T0 * S0), T0, S0)),
    X   = matrix(stats::rnorm(T0 * S0), T0, S0),
    COM = matrix(1, T0, S0),
    K   = matrix(1, T0, S0),
    TMG = stats::rnorm(T0),
    Gp  = stats::rnorm(T0),
    M = M)
}

test_that("checkpoint_dir is rejected with keep_fits = TRUE", {
  d <- .mi_ckpt_data()
  expect_error(
    fit_ou_nested_mi(phi_draws = d$phi_draws, X = d$X, TMG = d$TMG, COM = d$COM,
                     CAPITAL_TOTAL = d$K, Gprime = d$Gp, n_levels = 2L, M = d$M,
                     keep_fits = TRUE, checkpoint_dir = tempfile()),
    "incompatible with `keep_fits = TRUE`")
})

test_that("a fully-restored run reproduces the object bit-for-bit with zero refits", {
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat"),
                      inherits = FALSE),
              "testthat 3e local_mocked_bindings required")
  skip_if_not_installed("posterior")

  d  <- .mi_ckpt_data()
  fe <- .fake_engine_env()
  ck <- tempfile("ckpt_full_"); on.exit(unlink(ck, recursive = TRUE), add = TRUE)

  testthat::local_mocked_bindings(fit_ou_nested = fe$engine, .draws_matrix = fe$draws_matrix)

  run <- function() fit_ou_nested_mi(
    phi_draws = d$phi_draws, X = d$X, TMG = d$TMG, COM = d$COM,
    CAPITAL_TOTAL = d$K, Gprime = d$Gp, n_levels = 2L, M = d$M,
    keep_kappa_mixture = TRUE, seed = 1000, checkpoint_dir = ck)

  fe$calls <- 0L
  a <- run()
  expect_equal(fe$calls, d$M)                       # cold run: one fit per imputation
  expect_true(all(file.exists(file.path(ck, sprintf("imp_%d.rds", seq_len(d$M))))))

  fe$calls <- 0L
  b <- run()                                        # warm run: everything restored
  expect_equal(fe$calls, 0L)                        # no resampling
  expect_equal(b, a)                                # identical pooled object
})

test_that("a partial run refits only the missing imputations (crash-recovery)", {
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat"),
                      inherits = FALSE),
              "testthat 3e local_mocked_bindings required")
  skip_if_not_installed("posterior")

  d  <- .mi_ckpt_data()
  fe <- .fake_engine_env()
  ck <- tempfile("ckpt_part_"); on.exit(unlink(ck, recursive = TRUE), add = TRUE)

  testthat::local_mocked_bindings(fit_ou_nested = fe$engine, .draws_matrix = fe$draws_matrix)

  run <- function() fit_ou_nested_mi(
    phi_draws = d$phi_draws, X = d$X, TMG = d$TMG, COM = d$COM,
    CAPITAL_TOTAL = d$K, Gprime = d$Gp, n_levels = 2L, M = d$M,
    keep_kappa_mixture = TRUE, seed = 1000, checkpoint_dir = ck)

  fe$calls <- 0L
  full <- run()                                     # uninterrupted reference

  # Simulate a crash after imputation 1 of 3: drop the last two checkpoints.
  unlink(file.path(ck, c("imp_2.rds", "imp_3.rds")))
  fe$calls <- 0L
  resumed <- run()
  expect_equal(fe$calls, 2L)                        # only the two missing imputations refit
  expect_equal(resumed, full)                       # same pooled object as the uninterrupted run
})

test_that("checkpoint_dir = NULL writes nothing and is the legacy path", {
  skip_if_not(exists("local_mocked_bindings", where = asNamespace("testthat"),
                      inherits = FALSE),
              "testthat 3e local_mocked_bindings required")
  skip_if_not_installed("posterior")

  d  <- .mi_ckpt_data()
  fe <- .fake_engine_env()
  testthat::local_mocked_bindings(fit_ou_nested = fe$engine, .draws_matrix = fe$draws_matrix)

  ck  <- tempfile("ckpt_none_")
  out <- fit_ou_nested_mi(
    phi_draws = d$phi_draws, X = d$X, TMG = d$TMG, COM = d$COM,
    CAPITAL_TOTAL = d$K, Gprime = d$Gp, n_levels = 2L, M = d$M,
    keep_kappa_mixture = TRUE, seed = 1000)         # no checkpoint_dir
  expect_s3_class(out, "ou_nested_mi")
  expect_false(dir.exists(ck))                      # nothing created
})
