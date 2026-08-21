
test_that("Simulation with no covariates",{
  fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gengamma")
  sim <- simulate(fitg)
  expect_true(is.data.frame(sim))
})

test_that("Simulation with covariates",{
  ovarian$rx2 <- as.numeric(ovarian$rx == 2)
  fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ rx2, data = ovarian, dist="gamma")
  sim <- simulate(fitg, nsim=1)
  expect_true(is.data.frame(sim))
  sim <- simulate(fitg, nsim=2)
  expect_true(is.data.frame(sim))
  sim <- simulate(fitg, nsim=2, censtime=5)
  expect_true(is.data.frame(sim))
  sim <- simulate(fitg, nsim=2, censtime=seq(5, 10, length.out = nrow(ovarian)))
  expect_true(is.data.frame(sim))
})

test_that("Setting seed makes simulations reproducible",{
  sim1 <- simulate(fitg, nsim=2, censtime=5, seed=2)
  sim2 <- simulate(fitg, nsim=2, censtime=5, seed=2)
  expect_equal(sim1, sim2)
})

