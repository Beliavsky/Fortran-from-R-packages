test_that("rsdc_likelihood matches Hamilton (no X)", {
  y <- toy_residuals(50, 3); K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  P <- toy_P()

  # pack params: trans (2) + rho (N * K(K-1)/2 = 2*3=6)
  trans <- c(P[1,1], P[2,2])
  params <- c(trans, as.vector(t(rho)))
  nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)

  ref <- rsdc_hamilton(y = y, X = NULL, beta = NULL,
                       rho_matrix = rho, K = K, N = N, P = P)
  expect_equal(nll, -ref$log_likelihood, tolerance = 1e-10)
})

test_that("rsdc_likelihood matches Hamilton (TVTP)", {
  y <- toy_residuals(50, 3); K <- ncol(y); N <- 2
  rho <- toy_rho_matrix(K)
  X <- cbind(1, scale(seq_len(nrow(y))))
  p <- ncol(X)
  beta <- rbind(c(1.5, 0.1), c(0.8, -0.2))

  params <- c(as.vector(t(beta)), as.vector(t(rho)))
  nll <- rsdc_likelihood(params, y = y, exog = X, K = K, N = N)

  ref <- rsdc_hamilton(y = y, X = X, beta = beta,
                       rho_matrix = rho, K = K, N = N, P = NULL)
  expect_equal(nll, -ref$log_likelihood, tolerance = 1e-10)
})

test_that("rsdc_likelihood matches Hamilton (N=3 noX)", {
  set.seed(7)
  y <- toy_residuals(30, 2); K <- ncol(y); N <- 3
  rho_matrix <- matrix(c(0.1, 0.3, 0.5), nrow = N)   # 3x1 (K=2 => 1 pair)
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.8, 0.1,
                0.1, 0.2, 0.7), N, N, byrow = TRUE)
  # pack: [p_i1, p_i2] per row (p_i3 = 1 - p_i1 - p_i2), then rho row-wise
  trans_params <- c(P[1,1], P[1,2], P[2,1], P[2,2], P[3,1], P[3,2])
  params <- c(trans_params, as.vector(t(rho_matrix)))

  nll <- rsdc_likelihood(params, y = y, exog = NULL, K = K, N = N)
  ref <- rsdc_hamilton(y = y, X = NULL, beta = NULL,
                       rho_matrix = rho_matrix, K = K, N = N, P = P)
  expect_equal(nll, -ref$log_likelihood, tolerance = 1e-10)
})
