test_that("rsdc_hamilton works with fixed transition matrix (no X)", {
  y <- toy_residuals(T = 40, K = 3)
  K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  P <- toy_P()

  fit <- rsdc_hamilton(y = y, X = NULL, beta = NULL,
                       rho_matrix = rho, K = K, N = N, P = P)

  expect_type(fit$log_likelihood, "double")
  expect_true(is.finite(fit$log_likelihood))
  expect_equal(dim(fit$filtered_probs), c(N, nrow(y)))
  expect_equal(dim(fit$smoothed_probs), c(N, nrow(y)))

  # columns of probabilities sum to 1 and are in [0,1]
  cs_f <- colSums(fit$filtered_probs)
  cs_s <- colSums(fit$smoothed_probs)
  expect_true(all(abs(cs_f - 1) < 1e-8))
  expect_true(all(abs(cs_s - 1) < 1e-8))
  expect_true(all(fit$filtered_probs >= 0 & fit$filtered_probs <= 1))
  expect_true(all(fit$smoothed_probs  >= 0 & fit$smoothed_probs  <= 1))
})

test_that("rsdc_hamilton works with TVTP (X, beta)", {
  y <- toy_residuals(T = 40, K = 3)
  K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  # simple covariate: intercept + trend
  X <- cbind(1, scale(seq_len(nrow(y))))
  p <- ncol(X)
  # beta rows per regime; make regime 1 more persistent
  beta <- rbind(c(2,  0.0),
                c(1, -0.1))

  fit <- rsdc_hamilton(y = y, X = X, beta = beta,
                       rho_matrix = rho, K = K, N = N, P = NULL)
  expect_true(is.finite(fit$log_likelihood))
  expect_equal(dim(fit$filtered_probs), c(N, nrow(y)))
  expect_equal(dim(fit$smoothed_probs), c(N, nrow(y)))
})

test_that("rsdc_hamilton xi_init warm-start produces valid output", {
  set.seed(42)
  y <- toy_residuals(T = 30, K = 3)
  K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  P   <- toy_P()

  # Pass 1: IS slice — extract terminal filtered state
  res1    <- rsdc_hamilton(y = y[1:20, ], rho_matrix = rho, K = K, N = N, P = P)
  xi_term <- res1$filtered_probs[, 20]

  # Pass 2: OOS slice initialised from IS terminal state
  res2 <- rsdc_hamilton(y = y[21:30, ], rho_matrix = rho, K = K, N = N, P = P,
                        xi_init = xi_term)

  expect_equal(dim(res2$filtered_probs), c(N, 10))
  expect_equal(dim(res2$smoothed_probs), c(N, 10))
  expect_true(is.finite(res2$log_likelihood))
  expect_true(all(abs(colSums(res2$smoothed_probs) - 1) < 1e-10))
})
