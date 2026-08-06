test_that("rsdc_minvar produces valid weights and volatility", {
  skip_if_not_installed("quadprog")

  T <- 20; K <- 3
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  pairs <- toy_cor_pairs(K)

  # simple time-varying pairwise corr: small oscillations
  P <- length(pairs)
  pred_corr <- matrix(0, T, P)
  for (t in 1:T) pred_corr[t, ] <- c(0.2, 0.1, 0.05) + 0.01*sin(t/3)

  res <- rsdc_minvar(sigma_matrix = S,
                     value_cols = colnames(S),
                     predicted_corr = pred_corr,
                     y = y,
                     long_only = TRUE)

  expect_equal(dim(res$weights), c(T, K))
  # weights sum to ~1
  expect_true(all(abs(rowSums(res$weights) - 1) < 1e-6))
  # long-only honored
  expect_true(all(res$weights >= -1e-10))
  expect_true(is.finite(res$volatility))
})

test_that("rsdc_maxdiv produces valid outputs", {
  skip_if_not_installed("Rsolnp")

  T <- 20; K <- 3
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  pairs <- toy_cor_pairs(K); P <- length(pairs)
  pred_corr <- matrix(rep(c(0.15, 0.10, 0.05), each = T), ncol = P)

  res <- rsdc_maxdiv(sigma_matrix = S,
                     value_cols = colnames(S),
                     predicted_corr = pred_corr,
                     y = y,
                     long_only = TRUE)

  expect_equal(dim(res$weights), c(T, K))
  expect_equal(length(res$returns), T)
  expect_equal(length(res$diversification_ratios), T)
  expect_true(all(abs(rowSums(res$weights) - 1) < 1e-6))
  expect_true(all(res$weights >= -1e-10))
  expect_true(is.finite(res$mean_diversification))
  expect_true(is.finite(res$volatility))
})


test_that("rsdc_minvar reorders/warns when colnames(y) disagree with value_cols (L2)", {
  skip_if_not_installed("quadprog")
  T <- 15; K <- 3
  set.seed(9)
  S <- toy_sigma_matrix(T, K)                       # columns sigma1..sigmaK
  colnames(S) <- paste0("A", 1:K)
  pred_corr <- matrix(0.1, T, choose(K, 2))
  y <- matrix(rnorm(T * K, sd = 0.01), T, K)
  colnames(y) <- paste0("A", 1:K)

  # same names, shuffled order -> reordered with a warning; results must match
  # the correctly ordered call
  y_shuffled <- y[, c(2, 3, 1)]
  expect_warning(
    mv_shuf <- rsdc_minvar(S, paste0("A", 1:K), pred_corr, y_shuffled),
    "reordered")
  mv_ok <- rsdc_minvar(S, paste0("A", 1:K), pred_corr, y)
  expect_equal(mv_shuf$weights, mv_ok$weights)
  expect_equal(mv_shuf$volatility, mv_ok$volatility)

  # unrelated names -> assume positional order, but warn
  y_odd <- y; colnames(y_odd) <- paste0("Z", 1:K)
  expect_warning(rsdc_minvar(S, paste0("A", 1:K), pred_corr, y_odd),
                 "assumed")
})
