# The C++ log-likelihood used by rsdc_likelihood() must reproduce the pure-R
# Hamilton filter rsdc_hamilton() to numerical precision (added in 1.4-0).

test_that("rsdc_hamilton C++ engine matches the R engine (filtered/smoothed/loglik)", {
  set.seed(7)
  T <- 120; K <- 3
  y <- scale(matrix(rnorm(T * K), T, K))
  X <- cbind(1, as.numeric(scale(seq_len(T))))

  cmp <- function(N, mode) {
    rho <- matrix(seq(0.1, 0.7, length.out = N * K * (K - 1) / 2), nrow = N, byrow = TRUE)
    if (mode == "fixedP") {
      free <- 0.8 / (N - 1)
      P <- matrix(free, N, N); diag(P) <- 0.2 + 0.8 - free * (N - 1)
      P <- P / rowSums(P)
      a <- rsdc_hamilton(y, NULL, NULL, rho, K, N, P, engine = "r")
      b <- rsdc_hamilton(y, NULL, NULL, rho, K, N, P, engine = "cpp")
    } else if (mode == "uniform") {
      a <- rsdc_hamilton(y, NULL, NULL, rho, K, N, engine = "r")
      b <- rsdc_hamilton(y, NULL, NULL, rho, K, N, engine = "cpp")
    } else {
      beta <- matrix(rnorm(N * (if (N == 2) ncol(X) else (N - 1) * ncol(X)), sd = 0.2), nrow = N)
      a <- rsdc_hamilton(y, X, beta, rho, K, N, engine = "r")
      b <- rsdc_hamilton(y, X, beta, rho, K, N, engine = "cpp")
    }
    expect_equal(a$log_likelihood, b$log_likelihood, tolerance = 1e-8)
    expect_equal(a$filtered_probs, b$filtered_probs, tolerance = 1e-8)
    expect_equal(a$smoothed_probs, b$smoothed_probs, tolerance = 1e-8)
  }
  cmp(2, "fixedP"); cmp(3, "fixedP"); cmp(2, "uniform"); cmp(2, "tvtp"); cmp(3, "tvtp")
})

test_that("C++ log-likelihood matches the R Hamilton filter (noX and tvtp)", {
  set.seed(1)
  T <- 120; K <- 3; N <- 2
  y <- scale(matrix(rnorm(T * K), T, K))
  rho <- rbind(c(0.10, 0.05, 0.00), c(0.60, 0.40, 0.30))

  # Fixed transition (noX): params = (p11, p22, rho rows...)
  par_nox <- c(0.9, 0.8, as.vector(t(rho)))
  ll_cpp  <- -rsdc_likelihood(par_nox, y = y, exog = NULL, K = K, N = N)
  P       <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)
  ll_R    <- rsdc_hamilton(y, NULL, NULL, rho, K, N, P)$log_likelihood
  expect_equal(ll_cpp, ll_R, tolerance = 1e-8)

  # TVTP: params = (beta rows..., rho rows...)
  X <- cbind(1, as.numeric(scale(seq_len(T))))
  beta <- rbind(c(1.2, 0.0), c(0.8, -0.1))
  par_tv  <- c(as.vector(t(beta)), as.vector(t(rho)))
  ll_cpp2 <- -rsdc_likelihood(par_tv, y = y, exog = X, K = K, N = N)
  ll_R2   <- rsdc_hamilton(y, X, beta, rho, K, N)$log_likelihood
  expect_equal(ll_cpp2, ll_R2, tolerance = 1e-8)
})

test_that("C++ log-likelihood matches the R filter for N = 4 (noX and tvtp)", {
  set.seed(3)
  T <- 200; K <- 2; N <- 4
  y <- scale(matrix(rnorm(T * K), T, K))
  rho <- rbind(0.05, 0.35, 0.6, 0.85)                # N x C, C = 1

  # noX: each row a valid simplex; free = first N-1 = 3 entries
  Prows <- rbind(c(0.70, 0.10, 0.10, 0.10),
                 c(0.10, 0.70, 0.10, 0.10),
                 c(0.10, 0.10, 0.70, 0.10),
                 c(0.10, 0.10, 0.10, 0.70))
  trans <- as.vector(t(Prows[, 1:(N - 1), drop = FALSE]))   # row-major free entries
  par_nox <- c(trans, as.vector(t(rho)))
  ll_cpp <- -rsdc_likelihood(par_nox, y = y, exog = NULL, K = K, N = N)
  ll_R   <- rsdc_hamilton(y, NULL, NULL, rho, K, N, Prows)$log_likelihood
  expect_equal(ll_cpp, ll_R, tolerance = 1e-8)

  # tvtp: beta is N x (N-1)*p
  X <- cbind(1, as.numeric(scale(seq_len(T)))); p <- ncol(X)
  beta <- matrix(rnorm(N * (N - 1L) * p, sd = 0.2), nrow = N)
  par_tv <- c(as.vector(t(beta)), as.vector(t(rho)))
  ll_cpp2 <- -rsdc_likelihood(par_tv, y = y, exog = X, K = K, N = N)
  ll_R2   <- rsdc_hamilton(y, X, beta, rho, K, N)$log_likelihood
  expect_equal(ll_cpp2, ll_R2, tolerance = 1e-8)
})

test_that("C++ log-likelihood matches the R filter for N = 3", {
  set.seed(2)
  T <- 150; K <- 2; N <- 3
  y <- scale(matrix(rnorm(T * K), T, K))
  rho <- rbind(0.1, 0.5, 0.85)                       # N x C, C = 1
  X <- cbind(1, as.numeric(scale(seq_len(T)))); p <- ncol(X)
  beta <- rbind(c(0.5, 0, 0.2, 0), c(0.3, 0, 0.1, 0), c(0.4, 0, -0.1, 0))
  par_tv <- c(as.vector(t(beta)), as.vector(t(rho)))
  ll_cpp <- -rsdc_likelihood(par_tv, y = y, exog = X, K = K, N = N)
  ll_R   <- rsdc_hamilton(y, X, beta, rho, K, N)$log_likelihood
  expect_equal(ll_cpp, ll_R, tolerance = 1e-8)
})
