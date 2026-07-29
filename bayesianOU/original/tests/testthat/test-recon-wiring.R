# Input-validation and nesting-invariant tests for k_uncertainty = "recon".
# These do not sample; they exercise the R-side dispatch and the arithmetic that
# guarantees recon reduces to meas as sigma_K -> 0.

mk_data <- function(Tn = 24L, S = 2L) {
  set.seed(1)
  Gprime <- abs(stats::rnorm(Tn)) * 0.02 + 0.05
  k_cost <- matrix(stats::runif(Tn * S, 50, 150), Tn, S)
  K      <- matrix(stats::runif(Tn * S, 800, 1500), Tn, S)
  X      <- k_cost + K * matrix(Gprime, Tn, S)         # X = k + K G' (eq. 9)
  list(Y = matrix(stats::rnorm(Tn * S), Tn, S), X = X, k_cost = k_cost, K = K,
       Gprime = Gprime, TMG = stats::rnorm(Tn),
       COM = matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S), Tn = Tn, S = S)
}

test_that("recon requires n_levels >= 2", {
  # Since Level 3 (D-IMPL-10.1) recon is valid for n_levels >= 2 (the latent Phi
  # measurement block, widened from == 2 to >= 2); only single-level rejects it.
  d <- mk_data()
  expect_error(
    fit_ou_nested(d$Y, d$X, d$TMG, d$COM, d$K, n_levels = 1L,
                  k_uncertainty = "recon", k_cost = d$k_cost),
    "solo aplica con n_levels >= 2"
  )
})

test_that("recon requires k_cost", {
  d <- mk_data()
  # Force the backend check to fail fast: there is no Stan call because the
  # k_cost guard fires during data assembly (before .run_stan_ou).
  expect_error(
    suppressWarnings(fit_ou_nested(d$Y, d$X, d$TMG, d$COM, d$K, n_levels = 2L,
                                   Gprime = d$Gprime, k_uncertainty = "recon",
                                   k_cost = NULL)),
    "requiere `k_cost`"
  )
})

test_that("recon rejects non-positive CAPITAL_TOTAL", {
  d <- mk_data()
  Kbad <- d$K; Kbad[3, 1] <- 0
  expect_error(
    suppressWarnings(fit_ou_nested(d$Y, d$X, d$TMG, d$COM, Kbad, n_levels = 2L,
                                   Gprime = d$Gprime, k_uncertainty = "recon",
                                   k_cost = d$k_cost)),
    "must be strictly positive"
  )
})

test_that("nesting invariant: standardized reconstruction at K = K_hat equals zX", {
  # The Stan transformed-parameter Phi_anchor_eff at z_K = 0 is
  #   (k + K_hat G' - mu_X) / sd_X,
  # which must equal the meas anchor zX exactly, so recon nests meas.
  d <- mk_data()
  T_train <- floor(d$Tn * 0.70)
  zX <- zscore_train(d$X, T_train)
  recon_std <- sweep(sweep(d$k_cost + d$K * matrix(d$Gprime, d$Tn, d$S),
                           2, zX$mu), 2, zX$sd, "/")
  expect_equal(unname(recon_std), unname(zX$Mz), tolerance = 1e-12)
})
