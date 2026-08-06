test_that("B = 1 reproduces the historical single-draw behaviour", {
  set.seed(7)
  x <- matrix(rnorm(200 * 3), 200, 3)
  y <- matrix(rnorm(200 * 6), 200, 6)

  # the default must stay exactly what it was before B/seed existed
  a <- span_as(x, y)
  b <- span_as(x, y, control = list(B = 1L, seed = 123L))
  expect_identical(a, b)

  pv1 <- spantest:::f_getpv_batch(x, y, ks = 1/3, L = c(0, 2))
  pv2 <- spantest:::f_getpv_batch(x, y, ks = 1/3, L = c(0, 2), B = 1L, seed = 123L)
  expect_identical(pv1, pv2)
})

test_that("L = 0 is unaffected by B (there is no perturbation to merge)", {
  set.seed(8)
  x <- matrix(rnorm(150 * 2), 150, 2)
  y <- matrix(rnorm(150 * 4), 150, 4)

  p1 <- span_as(x, y, control = list(L = 0, B = 1L))
  p9 <- span_as(x, y, control = list(L = 0, B = 9L))
  expect_identical(p1, p9)
})

test_that("B > 1 merges draws and stabilises the L > 0 p-value", {
  set.seed(9)
  x <- matrix(rnorm(150 * 2), 150, 2)
  y <- matrix(rnorm(150 * 8), 150, 8)

  merged <- span_as(x, y, control = list(L = 2, B = 20L))
  expect_true(all(vapply(merged, function(p) is.finite(p) && p > 0 && p <= 1, logical(1))))

  # the merged value must not be pinned to any single draw
  draws <- vapply(1:20, function(s)
    span_as(x, y, control = list(L = 2, B = 1L, seed = s))$CCTad_L2_k1, numeric(1))
  expect_gt(stats::sd(draws), 0)                       # draws genuinely differ
  expect_gte(merged$CCTad_L2_k1, min(draws))           # merged is not the luckiest draw
})

test_that("seed is honoured and B is validated", {
  set.seed(10)
  x <- matrix(rnorm(120 * 2), 120, 2)
  y <- matrix(rnorm(120 * 3), 120, 3)

  expect_false(isTRUE(all.equal(span_as(x, y, control = list(L = 2, seed = 1L))$CCTad_L2_k1,
                                span_as(x, y, control = list(L = 2, seed = 999L))$CCTad_L2_k1)))
  expect_error(span_as(x, y, control = list(B = 0)))
  expect_error(span_as(x, y, control = list(B = NA)))
})

test_that("control validation rejects non-integer B and seed", {
  set.seed(12)
  x <- matrix(rnorm(120 * 2), 120, 2); y <- matrix(rnorm(120 * 3), 120, 3)
  expect_error(span_as(x, y, control = list(B = 1.9)))      # would silently truncate to 1
  expect_error(span_as(x, y, control = list(seed = 1.5)))
  expect_error(span_as(x, y, control = list(B = c(1, 2))))
  expect_silent(span_as(x, y, control = list(B = 3L, seed = -7L)))
})

test_that("degenerate inputs return NA rather than NaN or an error", {
  expect_true(is.na(spantest:::f_cauchypv(numeric(0))))     # every asset missing
  expect_false(is.nan(spantest:::f_cauchypv(numeric(0))))
  # singular benchmark: GL now matches the package's NA convention
  set.seed(13)
  x <- matrix(rnorm(100), 100, 1); x <- cbind(x, x)         # collinear
  y <- matrix(rnorm(300), 100, 3)
  expect_true(is.na(span_gl_a(x, y, control = list(totsim = 20))$pval_LMC))
  expect_true(is.na(span_gl_ad(x, y, control = list(totsim = 20))$pval_LMC))
  # too-short sample fails loudly instead of returning NaN
  expect_error(spantest:::f_ttest(matrix(rnorm(6), 6, 1), k = 1/3), "subseries")
})

test_that("span_grs and span_bj are the same test computed two ways", {
  # Britten-Jones (1999): regressing ones on raw returns reproduces the GRS
  # statistic exactly. Two independent code paths agreeing to machine precision
  # is a correctness check on both -- and the reason the two must never be
  # reported as if they were independent tests.
  set.seed(4242)
  worst <- 0
  for (i in 1:15) {
    T <- sample(120:300, 1); K <- sample(1:4, 1); N <- sample(1:6, 1)
    x <- matrix(rnorm(T * K), T, K); y <- matrix(rnorm(T * N), T, N)
    worst <- max(worst, abs(span_grs(x, y)$stat - span_bj(x, y)$stat))
  }
  expect_lt(worst, 1e-9)

  # span_f1 is close but genuinely different: it must NOT collapse onto them
  set.seed(7)
  x <- matrix(rnorm(250 * 3), 250, 3); y <- matrix(rnorm(250 * 5), 250, 5)
  expect_gt(abs(span_grs(x, y)$stat - span_f1(x, y)$stat), 1e-6)
})

test_that("the classical tests are exactly calibrated under a Gaussian iid null", {
  # A wrong degree of freedom or scaling factor shows up here and nowhere else.
  skip_on_cran()
  set.seed(20260802)
  p <- t(vapply(1:400, function(i) {
    R1 <- matrix(rnorm(250 * 3), 250, 3)
    B  <- matrix(runif(3 * 4), 3, 4); B <- sweep(B, 2, colSums(B), "/")
    R2 <- R1 %*% B + matrix(rnorm(250 * 4, sd = 0.02), 250, 4)
    c(span_grs(R1, R2)$pval, span_hk(R1, R2)$pval, span_f1(R1, R2)$pval,
      span_f2(R1, R2)$pval, span_km(R1, R2)$pval)
  }, numeric(5)))
  # generous band: 400 reps give a MC s.e. of about 1.1% at the 5% level
  expect_true(all(colMeans(p < 0.05) > 0.02 & colMeans(p < 0.05) < 0.09))
})

test_that("each test responds to the hypothesis it claims and ignores the other", {
  # A test of the WRONG null is still correctly sized when both restrictions
  # hold, so calibration alone cannot catch it. Only discrimination can: drive
  # delta away from zero with alpha fixed at zero, and vice versa.
  skip_on_cran()
  set.seed(20260803)
  T <- 250L; K <- 3L; N <- 4L; NREP <- 300L
  gen <- function(alpha, colsum) {
    R1 <- matrix(rnorm(T * K), T, K)
    B  <- matrix(runif(K * N), K, N); B <- sweep(B, 2, colSums(B), "/") * colsum
    list(R1 = R1, R2 = sweep(R1 %*% B, 2, alpha, "+") + matrix(rnorm(T * N, sd = 0.05), T, N))
  }
  rate <- function(f, alpha, colsum)
    mean(vapply(seq_len(NREP), function(i) { d <- gen(alpha, colsum); f(d$R1, d$R2) },
                numeric(1)) < 0.05)

  km <- function(a, b) span_km(a, b)$pval
  f2 <- function(a, b) span_f2(a, b)$pval
  gr <- function(a, b) span_grs(a, b)$pval

  # delta tests must fire on delta and stay quiet on alpha
  expect_gt(rate(km, rep(0, N),    1.25), 0.90)   # delta != 0  -> reject
  expect_lt(rate(km, rep(0.03, N), 1.00), 0.15)   # alpha != 0  -> do not
  expect_gt(rate(f2, rep(0, N),    1.25), 0.90)
  expect_lt(rate(f2, rep(0.03, N), 1.00), 0.15)
  # and the alpha test the other way round
  expect_gt(rate(gr, rep(0.03, N), 1.00), 0.90)
  expect_lt(rate(gr, rep(0, N),    1.25), 0.15)
})
