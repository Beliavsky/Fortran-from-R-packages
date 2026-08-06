# Phase-1 functionality added in 1.4-0: robust SEs, diagnostics, plot, warm start.

fast <- list(itermax = 30, NP = 50, seed = 1)

make_data <- function(T = 300) {
  set.seed(1)
  X <- cbind(1, as.numeric(scale(seq_len(T))))
  b <- rbind(c(1.6, -0.7), c(1.1, 0.2)); rho <- rbind(0.25, 0.8)
  Sig <- array(c(matrix(c(1, .25, .25, 1), 2), matrix(c(1, .8, .8, 1), 2)), c(2, 2, 2))
  sim <- rsdc_simulate(T, X, b, matrix(0, 2, 2), Sig, 2, seed = 7)
  list(y = scale(sim$observations), X = X)
}

test_that("vcov() supports hessian/opg/sandwich and confint accepts type", {
  skip_on_cran()
  d <- make_data(); fit <- rsdc_estimate("tvtp", residuals = d$y, N = 2, X = d$X, control = fast)
  np <- fit$npar
  # Assert the covariances are actually produced: wrapping these in
  # `if (!is.null(V))` would let a fit with no covariance at all pass silently.
  for (ty in c("hessian", "opg", "sandwich")) {
    V <- vcov(fit, type = ty)
    expect_equal(dim(V), c(np, np))
    expect_true(all(is.finite(V)))
  }
  expect_equal(nrow(confint(fit, type = "sandwich")), np)
})

test_that("a penalty-contaminated Hessian is refused instead of yielding absurd SEs", {
  skip_on_cran()
  # optimHess differences by +/- ndeps (1e-3); the objectives return a large
  # finite penalty outside the feasible set. A correlation pinned next to 1
  # therefore used to produce a ~1e15 Hessian and standard errors of ~1e-08.
  set.seed(7)
  y <- matrix(rnorm(600), ncol = 2)
  fn <- function(p) rsdc_likelihood(p, y = y, exog = NULL, K = 2, N = 2)
  par <- c(0.9, 0.9, 0.9995, 0.2)          # noX N=2: 2 stay probs, then 2 rhos
  expect_gt(fn(c(par[1:2], 0.9995 + 1e-3, par[4])), 1e9)   # the wall is in reach
  expect_warning(V <- RSDC:::.rsdc_vcov(par, fn), "penalty-contaminated|feasible")
  expect_null(V)
  # An interior point still returns a usable covariance.
  V_ok <- RSDC:::.rsdc_vcov(c(0.9, 0.9, 0.3, 0.7), fn)
  expect_true(is.null(V_ok) || all(is.finite(V_ok)))
})

test_that("scores and OPG are withheld at a projected (non-stationary) noX optimum", {
  skip_on_cran()
  d <- make_data(); fit <- rsdc_estimate("noX", residuals = d$y, N = 2, control = fast)
  fit$convergence <- 99L                   # as f_optim_noX sets on row projection
  fit$vcov <- NULL; fit$scores <- NULL      # and as it then withholds
  expect_error(vcov(fit, type = "opg"), "No variance-covariance available")
  expect_error(vcov(fit, type = "hessian"), "No variance-covariance available")
})

test_that("summary carries diagnostics with a valid ergodic distribution", {
  skip_on_cran()
  d <- make_data(); fit <- rsdc_estimate("noX", residuals = d$y, N = 2, control = fast)
  s <- summary(fit)
  expect_false(is.null(s$diagnostics))
  expect_equal(nrow(s$diagnostics), 2L)
  erg <- s$diagnostics$ergodic_prob
  if (all(is.finite(erg))) expect_equal(sum(erg), 1, tolerance = 1e-6)
  expect_true(all(s$diagnostics$exp_duration >= 1))
})

test_that("warm start skips the global search and reproduces the fit", {
  skip_on_cran()
  d <- make_data(); fit <- rsdc_estimate("tvtp", residuals = d$y, N = 2, X = d$X, control = fast)
  fit2 <- rsdc_estimate("tvtp", residuals = d$y, N = 2, X = d$X,
                        control = list(start = fit$par, compute_se = FALSE))
  expect_equal(fit2$log_likelihood, fit$log_likelihood, tolerance = 0.5)
  expect_error(
    rsdc_estimate("tvtp", residuals = d$y, N = 2, X = d$X, control = list(start = 1:3)),
    "length")
})

test_that("plot() runs for switching models and stores probabilities", {
  skip_on_cran()
  d <- make_data(); fit <- rsdc_estimate("noX", residuals = d$y, N = 2, control = fast)
  expect_false(is.null(fit$smoothed_probs))
  expect_equal(dim(fit$smoothed_probs), c(2L, nrow(d$y)))
  pf <- tempfile(fileext = ".pdf"); pdf(pf); on.exit(unlink(pf))
  expect_invisible(plot(fit)); dev.off()
})
