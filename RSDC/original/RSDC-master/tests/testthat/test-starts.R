# rsdc_starts() and the rsdc_starts branch of rsdc_estimate() (added in 1.6-0).

sim_data <- function(n = 300, K = 3, seed = 42) {
  rho <- rbind(seq(0.05, by = 0.05, length.out = K * (K - 1) / 2),
               seq(0.60, by = 0.04, length.out = K * (K - 1) / 2))
  Sig <- array(0, dim = c(K, K, 2))
  for (m in 1:2) {
    R <- diag(K); R[lower.tri(R)] <- rho[m, ]; R[upper.tri(R)] <- t(R)[upper.tri(R)]
    Sig[, , m] <- R
  }
  rsdc_simulate(n = n, X = matrix(1, n, 1), beta = matrix(qlogis(0.9), 2, 1),
                mu = matrix(0, 2, K), sigma = Sig, N = 2, seed = seed)
}

test_that("rsdc_starts returns feasible, correctly sized starts", {
  sim <- sim_data()
  y <- sim$observations
  X <- cbind(1, as.numeric(scale(seq_len(nrow(y)))))

  for (spec in list(list(m = "noX",  N = 2), list(m = "noX",  N = 3),
                    list(m = "tvtp", N = 2), list(m = "tvtp", N = 3),
                    list(m = "const", N = 1))) {
    st <- rsdc_starts(y, N = spec$N, method = spec$m,
                      X = if (spec$m == "tvtp") X else NULL,
                      window = 60, n_starts = 3)
    expect_s3_class(st, "rsdc_starts")
    expect_equal(nrow(st$starts), 3L)
    expect_true(all(is.finite(st$loglik0)))
    # start length matches the model's parameter count
    fit1 <- suppressWarnings(rsdc_estimate(
      spec$m, residuals = y, N = spec$N,
      X = if (spec$m == "tvtp") X else NULL,
      control = list(start = st$starts[1, ], compute_se = FALSE)))
    expect_equal(ncol(st$starts), fit1$npar)
    expect_identical(colnames(st$starts), names(fit1$coefficients))
  }
})

test_that("shrink overrides n_starts and both are validated", {
  y <- sim_data()$observations
  st <- rsdc_starts(y, N = 2, method = "noX", window = 60,
                    shrink = c(1, 0.5))
  expect_equal(nrow(st$starts), 2L)
  expect_equal(st$shrink, c(1, 0.5))
  expect_error(rsdc_starts(y, N = 2, method = "noX", window = 60, shrink = c(0, 1)),
               "shrink")
  expect_error(rsdc_starts(y, N = 2, method = "noX", window = 60, n_starts = 0),
               "n_starts")
  expect_error(rsdc_starts(y[1:50, ], N = 2, method = "noX", window = 60),
               "too short")
})

test_that("rsdc_estimate accepts an rsdc_starts object and keeps the best fit", {
  y <- sim_data()$observations
  st <- rsdc_starts(y, N = 2, method = "noX", window = 60, n_starts = 2)

  fit <- suppressWarnings(rsdc_estimate("noX", residuals = y, N = 2,
                                        control = list(start = st)))
  expect_s3_class(fit, "rsdc_fit")
  expect_length(fit$start_logliks, 2L)
  expect_equal(dim(fit$start_pars), c(2L, fit$npar))
  expect_identical(colnames(fit$start_pars), names(fit$coefficients))

  # level 2 reproduces a manual level-1 loop over the same starts
  manual <- lapply(seq_len(nrow(st$starts)), function(i)
    suppressWarnings(rsdc_estimate("noX", residuals = y, N = 2,
                                   control = list(start = st$starts[i, ]))))
  lls <- vapply(manual, function(f) f$log_likelihood, numeric(1))
  expect_equal(sort(fit$start_logliks), sort(lls), tolerance = 1e-6)
  expect_equal(fit$log_likelihood, max(lls), tolerance = 1e-6)
})

test_that("incompatible rsdc_starts objects are rejected with clear errors", {
  sim <- sim_data()
  y <- sim$observations
  X <- cbind(1, as.numeric(scale(seq_len(nrow(y)))))
  st <- rsdc_starts(y, N = 2, method = "noX", window = 60, n_starts = 2)

  expect_error(rsdc_estimate("tvtp", residuals = y, N = 2, X = X,
                             control = list(start = st)), "method")
  expect_error(rsdc_estimate("noX", residuals = y, N = 3,
                             control = list(start = st)), "N = 2")
  expect_error(rsdc_estimate("noX", residuals = y[, 1:2], N = 2,
                             control = list(start = st)), "K = 3")
})

test_that("cores > 1 reproduces the sequential result exactly", {
  y <- sim_data()$observations
  st <- rsdc_starts(y, N = 2, method = "noX", window = 60, n_starts = 2)

  f1 <- suppressWarnings(rsdc_estimate("noX", residuals = y, N = 2,
                                       control = list(start = st, cores = 1)))
  f2 <- suppressWarnings(rsdc_estimate("noX", residuals = y, N = 2,
                                       control = list(start = st, cores = 2)))
  expect_equal(f1$par, f2$par)
  expect_equal(f1$start_logliks, f2$start_logliks)
})
