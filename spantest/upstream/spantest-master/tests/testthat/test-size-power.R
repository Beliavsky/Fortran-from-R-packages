# Size (under H0) and power (under H1) of the spanning tests, using the
# package's own span_simulate() as the data-generating process. Heavy, so
# skipped on CRAN.

test_that("spanning tests have approximately nominal size under H0", {
  skip_on_cran()
  reps <- 150L
  pmat <- vapply(seq_len(reps), function(s) {
    set.seed(s)
    sim <- span_simulate(n = 250, K = 3, N = 6, ncp = 0, dgp = 1)
    x <- sim$R1; y <- sim$R2
    c(GRS = span_grs(x, y)$pval,
      HK  = span_hk(x, y)$pval,
      F1  = span_f1(x, y)$pval,
      F2  = span_f2(x, y)$pval,
      KM  = span_km(x, y)$pval,
      PY  = span_py(x, y)$pval,
      GLa = unname(span_gl_a(x, y, control = list(totsim = 100))$pval_LMC))
  }, numeric(7))
  size <- rowMeans(pmat < 0.05, na.rm = TRUE)
  # ~150 reps => MC se ~ 0.018; require size within 0.07 of nominal (catches
  # gross distortions like ~0 or ~0.5 without being flaky).
  expect_true(all(abs(size - 0.05) < 0.07),
              info = paste(names(size), round(size, 3), collapse = "; "))
})

test_that("spanning tests gain power under a mean-variance alternative", {
  skip_on_cran()
  reps <- 150L
  pmat <- vapply(seq_len(reps), function(s) {
    set.seed(1000L + s)
    sim <- span_simulate(n = 250, K = 3, N = 6, ncp = 0.4, dgp = 1)
    x <- sim$R1; y <- sim$R2
    c(GRS = span_grs(x, y)$pval,
      HK  = span_hk(x, y)$pval,
      F1  = span_f1(x, y)$pval,
      KM  = span_km(x, y)$pval)
  }, numeric(4))
  power <- rowMeans(pmat < 0.05, na.rm = TRUE)
  expect_true(all(power > 0.5),
              info = paste(names(power), round(power, 3), collapse = "; "))
})

test_that("F-based tests return NA on singular (collinear) benchmarks", {
  set.seed(1)
  X  <- matrix(rnorm(250 * 2), 250, 2)
  R1 <- cbind(X, X[, 1])                 # perfectly collinear benchmark column
  R2 <- matrix(rnorm(250 * 4), 250, 4)
  expect_true(is.na(span_grs(R1, R2)$pval))
  expect_true(is.na(span_f2(R1, R2)$pval))
  expect_true(is.na(span_py(R1, R2)$pval))
})
