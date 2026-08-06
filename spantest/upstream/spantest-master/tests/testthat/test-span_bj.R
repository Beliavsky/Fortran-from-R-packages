test_that("span_bj returns expected structure and values", {
  set.seed(321)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  result <- span_bj(R1, R2)

  expect_type(result, "list")
  expect_named(result, c("pval", "stat", "H0"))
  expect_true(is.numeric(result$stat) && length(result$stat) == 1)
  expect_true(is.numeric(result$pval) && result$pval >= 0 && result$pval <= 1)
  expect_equal(result$H0, "alpha = 0")
})

test_that("span_bj has power against an alpha alternative", {
  # Test assets carry a large intercept (alpha = 1): the tangency portfolio is
  # not spanned by the benchmarks, so the test must reject.
  set.seed(321)
  T <- 250; K <- 3; N <- 3
  R1 <- matrix(rnorm(T * K), T, K)
  B  <- matrix(runif(K * N, 0.2, 1), K, N)
  B  <- sweep(B, 2, colSums(B), "/")
  R2 <- 1 + R1 %*% B + matrix(rnorm(T * N, sd = 0.3), T, N)

  expect_lt(span_bj(R1, R2)$pval, 0.01)
})

test_that("span_bj matches the canonical Britten-Jones F statistic", {
  set.seed(1)
  R1 <- matrix(rnorm(200 * 3), 200, 3)
  R2 <- matrix(rnorm(200 * 2), 200, 2)

  R <- cbind(R1, R2); Tn <- nrow(R); K <- ncol(R1); N <- ncol(R2); ones <- rep(1, Tn)
  rssU <- sum(lm.fit(R,  ones)$residuals^2)
  rssR <- sum(lm.fit(R1, ones)$residuals^2)
  stat_ref <- ((rssR - rssU) / N) / (rssU / (Tn - K - N))

  expect_equal(span_bj(R1, R2)$stat, stat_ref, tolerance = 1e-8)
})
