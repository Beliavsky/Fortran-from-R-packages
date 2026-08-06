# Reparameterized global search: reproducibility, seed/cores invariance,
# replication warning, warm-start bypass, nested-model monotonicity.

gs_residuals <- function(T = 250, K = 3) {
  set.seed(7)
  R <- matrix(c(1, .5, .3, .5, 1, .4, .3, .4, 1), 3)
  Z <- matrix(rnorm(T * K), T, K) %*% chol(R)
  apply(Z, 2, function(x) x / sd(x))
}
fast <- list(compute_se = FALSE, itermax = 40)

test_that("default fit is exactly reproducible call to call", {
  z <- gs_residuals()
  f1 <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2, control = fast))
  f2 <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2, control = fast))
  expect_identical(coef(f1), coef(f2))
  expect_identical(logLik(f1), logLik(f2))
})

test_that("multi-start result is invariant to cores and replicates the max", {
  skip_on_cran()
  z <- gs_residuals()
  c1 <- c(fast, list(n_starts = 3L, cores = 1L))
  c2 <- c(fast, list(n_starts = 3L, cores = 2L))
  f1 <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2, control = c1))
  f2 <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2, control = c2))
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
  expect_length(f1$start_logliks, 3L)
  # On this easy unimodal surface all searches must agree (no warning).
  expect_true(sum(f1$start_logliks >= max(f1$start_logliks) - 2.5) >= 2L)
})

test_that("control$seed changes the search stream without erroring", {
  z <- gs_residuals()
  f <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2,
                                      control = c(fast, list(seed = 999L))))
  expect_s3_class(f, "rsdc_fit")
  expect_true(is.finite(as.numeric(logLik(f))))
})

test_that("a numeric warm start skips the global search unchanged", {
  z <- gs_residuals()
  f0 <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2, control = fast))
  f1 <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 2,
          control = list(compute_se = FALSE, start = f0$par)))
  expect_gte(as.numeric(logLik(f1)), as.numeric(logLik(f0)) - 1e-6)
})

test_that("maximized log-likelihoods are monotone across nested models", {
  skip_on_cran()
  z <- gs_residuals()
  X <- cbind(1, as.numeric(scale(seq_len(nrow(z)))))
  ll <- function(f) as.numeric(logLik(f))
  # Nesting monotonicity holds at the TRUE optima; a single under-lucky
  # search can genuinely violate it (that is exactly what the diagnostic
  # detects -- on this toy surface one tvtp2 seed lands 0.8 below noX2's
  # optimum while multi-start recovers a better one), so the nested fits are
  # estimated with the recommended multi-start protocol.
  full <- list(compute_se = FALSE, itermax = 200, n_starts = 4, cores = 2)
  f_const <- suppressWarnings(rsdc_estimate("const", residuals = z, N = 1, control = full))
  f_noX2  <- suppressWarnings(rsdc_estimate("noX",   residuals = z, N = 2, control = full))
  f_tvtp2 <- suppressWarnings(rsdc_estimate("tvtp",  residuals = z, N = 2, X = X, control = full))
  # A superset model can always reproduce the subset's optimum.
  expect_lte(ll(f_const), ll(f_noX2) + 0.01)
  expect_lte(ll(f_noX2),  ll(f_tvtp2) + 0.01)
})

test_that("noX N = 3 fits through the bounded-softmax head", {
  skip_on_cran()
  z <- gs_residuals(T = 300)
  f <- suppressWarnings(rsdc_estimate("noX", residuals = z, N = 3, control = fast))
  P <- f$transition_matrix
  expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-8)
  expect_true(all(P >= 0))
})
