# Features added in 1.5-0: multi-start, Viterbi, multi-step forecast,
# correlation bands, broom tidiers, and the ggplot2 autoplot method.

test_that("multi-start keeps the best fit and records the start log-likelihoods", {
  skip_on_cran()
  set.seed(1); y <- scale(matrix(rnorm(300 * 2), 300, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2,
                     control = list(n_starts = 3L, itermax = 40))
  expect_s3_class(f, "rsdc_fit")
  expect_length(f$start_logliks, 3L)
  expect_equal(f$log_likelihood, max(f$start_logliks), tolerance = 1e-6)
})

test_that("rsdc_viterbi returns a valid most-likely regime path", {
  skip_on_cran()
  set.seed(2); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  v <- rsdc_viterbi(f)
  expect_length(v, 200L)
  expect_true(all(v %in% 1:2))
})

test_that("rsdc_viterbi recovers a well-separated regime path", {
  skip_on_cran()
  set.seed(11)
  TT <- 400L
  X <- cbind(1, rep(0, TT))                       # constant covariate -> constant P
  lg <- function(p) log(p / (1 - p))
  beta <- rbind(c(lg(0.97), 0), c(lg(0.97), 0))   # persistent stay probs
  Sig  <- array(c(matrix(c(1, 0.0, 0.0, 1), 2),   # regime 1: low correlation
                  matrix(c(1, 0.95, 0.95, 1), 2)), # regime 2: high correlation
                c(2, 2, 2))
  s <- rsdc_simulate(TT, X, beta, matrix(0, 2, 2), Sig, N = 2, seed = 3)
  f <- rsdc_estimate("tvtp", residuals = scale(s$observations), N = 2, X = X,
                     control = list(itermax = 40, seed = 1))
  v <- rsdc_viterbi(f)
  expect_true(mean(v == s$states) > 0.75)         # MAP path tracks the true states
})

test_that("rsdc_forecast_ahead regime probabilities converge to the ergodic distribution (noX)", {
  skip_on_cran()
  set.seed(12); y <- scale(matrix(rnorm(300 * 2), 300, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  fa <- rsdc_forecast_ahead(f, horizon = 300L)
  P  <- f$transition_matrix
  ev <- eigen(t(P)); pic <- Re(ev$vectors[, which.min(abs(ev$values - 1))]); pic <- pic / sum(pic)
  expect_equal(as.numeric(fa$regime_probs[300, ]), pic, tolerance = 1e-3)
})

test_that("rsdc_forecast_ahead returns h-step regime and correlation forecasts", {
  skip_on_cran()
  set.seed(3); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  fa <- rsdc_forecast_ahead(f, horizon = 5L)
  expect_equal(dim(fa$regime_probs), c(5L, 2L))
  expect_equal(dim(fa$predicted_correlations), c(5L, 1L))
  expect_true(all(abs(rowSums(fa$regime_probs) - 1) < 1e-8))
})

test_that("rsdc_corr_bands returns ordered pointwise bands", {
  skip_on_cran()
  set.seed(4); y <- scale(matrix(rnorm(250 * 2), 250, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2,
                     control = list(itermax = 40, compute_se = TRUE))
  b <- rsdc_corr_bands(f, B = 50L, seed = 1)
  expect_length(b, 1L)                       # one asset pair for K = 2
  m <- b[[1]]
  expect_equal(colnames(m), c("fit", "lower", "upper"))
  expect_equal(nrow(m), nrow(y))
  expect_true(all(m[, "lower"] <= m[, "upper"]))
})

test_that("broom tidiers (tidy/glance/augment) work", {
  skip_on_cran()
  set.seed(5); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  td <- generics::tidy(f)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in% names(td)))
  gl <- generics::glance(f)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("logLik", "AIC", "BIC", "nobs") %in% names(gl)))
  au <- generics::augment(f)
  expect_equal(nrow(au), 200L)
  expect_true(".state" %in% names(au))
})

test_that("vcov/confint accept the simulation-based type = 'bootstrap'", {
  skip_on_cran()
  set.seed(7); y <- scale(matrix(rnorm(250 * 2), 250, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  np <- length(coef(f))
  Vb <- vcov(f, type = "bootstrap", B = 30, seed = 1)
  expect_equal(dim(Vb), c(np, np))
  expect_true(all(diag(Vb) >= 0))
  ci <- confint(f, type = "bootstrap", B = 30, seed = 1)
  expect_equal(nrow(ci), np)
  expect_true(all(ci[, 1] <= ci[, 2]))
})

test_that("autoplot returns a ggplot for switching models", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  set.seed(6); y <- scale(matrix(rnorm(200 * 2), 200, 2))
  f <- rsdc_estimate("noX", residuals = y, N = 2, control = list(itermax = 40))
  p <- ggplot2::autoplot(f)
  expect_s3_class(p, "ggplot")
})

test_that("rsdc_forecast_ahead validates the width of X_future (L3)", {
  skip_on_cran()
  set.seed(5); Tn <- 120
  y <- scale(matrix(rnorm(Tn * 2), Tn, 2))
  X <- cbind(1, as.numeric(scale(seq_len(Tn))))
  f <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X,
                     control = list(itermax = 40, compute_se = FALSE))
  # wrong number of columns must error, not silently mis-multiply
  expect_error(rsdc_forecast_ahead(f, horizon = 3L, X_future = matrix(1, 3, 3)),
               "X_future must have 2 column")
  # correct width works
  fa <- rsdc_forecast_ahead(f, horizon = 3L, X_future = cbind(1, c(0, 0.1, 0.2)))
  expect_equal(dim(fa$regime_probs), c(3L, 2L))
})
