# Correctness guards for the net-new geometry engine (R/geometry.R). Pure R
# (no Stan backend), so always-on. They check that (i) the Euclidean HMC samples
# a standard normal, (ii) the SoftAbs Riemannian metric equals the precision on
# a Gaussian (analytic Hessian) and yields a near-exact preconditioner
# (acceptance ~1, recovered variance), and (iii) the metrics are SPD. The
# Riemannian sampler must stay correct regardless of the metric because the
# metric is a preconditioner, not part of the target (Metropolis exact).

test_that("Euclidean HMC samples a standard normal", {
  tgt <- ou_geom_target(log_prob = function(th) -0.5 * sum(th^2),
                        grad_log_prob = function(th) -th, dim = 3L)
  fit <- ou_geom_hmc(tgt, epsilon = 0.3, L = 12L, n_iter = 800L,
                     n_warmup = 400L, seed = 1)
  expect_gt(fit$accept_rate, 0.7)
  expect_identical(fit$n_divergent, 0L)
  expect_lt(abs(mean(fit$draws[, 1])), 0.25)
  expect_lt(abs(sd(fit$draws[, 1]) - 1), 0.3)
})

test_that("SoftAbs Riemannian metric equals the precision on a Gaussian", {
  vars <- c(1, 100); iv <- 1 / vars
  tgt <- ou_geom_target(log_prob = function(th) -0.5 * sum(th^2 * iv),
                        grad_log_prob = function(th) -th * iv,
                        hessian = function(th) -diag(iv), dim = 2L)
  m <- ou_geom_metric_riemannian(tgt, curvature = "softabs")
  M <- m$mass(c(0, 0))
  # SoftAbs of the Hessian of -log pi (= diag(iv)) is ~ diag(iv) = the precision.
  expect_equal(diag(M), iv, tolerance = 1e-3)
  expect_true(all(eigen(M, symmetric = TRUE)$values > 0))   # SPD
  fit <- ou_geom_hmc(tgt, metric = m, epsilon = 0.2, L = 10L,
                     n_iter = 400L, n_warmup = 200L, seed = 3)
  expect_gt(fit$accept_rate, 0.9)
  expect_identical(fit$n_divergent, 0L)
  expect_lt(abs(sd(fit$draws[, 2]) - sqrt(vars[2])), 3)     # recovers wide sd
})

test_that("Euclidean dense metric is validated and the identity is the default", {
  m0 <- ou_geom_metric_euclidean(dim = 4L)
  expect_false(m0$position_dependent)
  expect_equal(m0$mass(rep(0, 4)), diag(4))
  expect_error(ou_geom_metric_euclidean(dim = 2L, M = matrix(c(1, 2, 2, 1), 2)),
               "positive-definite")
  md <- ou_geom_metric_euclidean(M = c(1, 4, 9))
  expect_equal(diag(md$inv_mass(rep(0, 3))), 1 / c(1, 4, 9))
})
