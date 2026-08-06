# Correctness of internal helpers and compiled kernels -----------------------

test_that("f_getpv_batch equals a column loop of f_getpv", {
  # This equivalence is what the vectorised span_as() speedup relies on.
  set.seed(1)
  Tn <- 200L; K <- 3L; N <- 6L
  x <- matrix(rnorm(Tn * K), Tn, K)
  Y <- matrix(rnorm(Tn * N), Tn, N)
  ks <- c(1 / 3); L <- c(0, 2)

  batch <- spantest:::f_getpv_batch(x, Y, ks = ks, L = L)
  loop  <- lapply(seq_len(N), function(j) spantest:::f_getpv(Y[, j], x, ks = ks, L = L))

  for (nm in names(batch)) {
    col_vals <- vapply(loop, function(z) unname(z[nm]), numeric(1))
    expect_equal(unname(batch[[nm]]), col_vals, info = nm)
  }
})

test_that("garch_filter (C++) matches the plain R GARCH(1,1) recursion", {
  gR <- function(z, om, al, be) {
    n <- length(z); eps <- numeric(n); s2 <- om / (1 - al - be)
    for (t in seq_len(n)) { e <- sqrt(s2) * z[t]; eps[t] <- e; s2 <- om + al * e * e + be * s2 }
    eps
  }
  set.seed(1); z <- rnorm(500)
  expect_equal(spantest:::garch_filter(z, 0.1, 0.1, 0.8), gR(z, 0.1, 0.1, 0.8))
  expect_equal(spantest:::garch_filter(z, 0.05, 0.1, 0.85), gR(z, 0.05, 0.1, 0.85))
})

test_that("f_cauchypv returns the common value when all p-values are equal", {
  expect_equal(spantest:::f_cauchypv(rep(0.5, 5)), 0.5)
  expect_equal(spantest:::f_cauchypv(rep(0.2, 4)), 0.2, tolerance = 1e-8)
  expect_equal(spantest:::f_cauchypv(rep(0.8, 3)), 0.8, tolerance = 1e-8)
})

test_that("f_rsstd is zero-mean, unit-variance, and skews with xi", {
  set.seed(1)
  x <- spantest:::f_rsstd(2e5, nu = 6, xi = 0.9)   # xi < 1 => left-skewed
  expect_lt(abs(mean(x)), 0.02)
  expect_lt(abs(var(x) - 1), 0.05)
  skew <- function(u) mean((u - mean(u))^3) / mean((u - mean(u))^2)^1.5
  expect_lt(skew(x), 0)
  set.seed(2)
  s <- spantest:::f_rsstd(2e5, nu = 6, xi = 1)      # symmetric
  expect_lt(abs(skew(s)), 0.05)
})

test_that("f_ranklex ranks the last element with the documented convention", {
  uu <- c(0.1, 0.2, 0.3)
  expect_identical(spantest:::f_ranklex(c(1, 2, 3), uu), 3)  # last is largest
  expect_identical(spantest:::f_ranklex(c(3, 2, 1), uu), 1)  # last is smallest
  expect_true(is.na(spantest:::f_ranklex(5, 0.5)))           # length 1 -> NA
  expect_true(is.na(spantest:::f_ranklex(numeric(0), numeric(0))))
})
