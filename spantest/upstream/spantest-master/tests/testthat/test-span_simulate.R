test_that("span_simulate returns correct dimensions", {
  set.seed(1)
  sim <- span_simulate(n = 200, K = 3, N = 8, ncp = 0)
  expect_type(sim, "list")
  expect_equal(dim(sim$R1), c(200L, 3L))
  expect_equal(dim(sim$R2), c(200L, 8L))
  expect_true(all(is.finite(sim$R1)))
  expect_true(all(is.finite(sim$R2)))
})

test_that("span_simulate is reproducible given the RNG seed", {
  set.seed(42); a <- span_simulate(120, 2, 5, ncp = 0.1, dgp = 4)
  set.seed(42); b <- span_simulate(120, 2, 5, ncp = 0.1, dgp = 4)
  expect_identical(a, b)
})

test_that("all DGP presets 1:12 run and return finite returns", {
  for (d in 1:12) {
    set.seed(d)
    sim <- span_simulate(150, 3, 6, ncp = 0, dgp = d)
    expect_equal(dim(sim$R2), c(150L, 6L), info = paste("dgp", d))
    expect_true(all(is.finite(sim$R2)), info = paste("dgp", d))
  }
})

test_that("dgp presets encode the documented raw vs standardised t variance", {
  # DGP 2 (raw t_5): variance ~ 5/3; DGP 5 (standardised t_4 GARCH): variance ~ 1
  set.seed(1); v_raw <- var(as.numeric(span_simulate(1e5, 1, 1, dgp = 2)$R1[, 1]))
  set.seed(1); v_std <- var(as.numeric(span_simulate(1e5, 1, 1, dgp = 5)$R1[, 1]))
  expect_gt(v_raw, 1.4)          # raw t_5 has variance 5/3 ~ 1.67
  expect_lt(abs(v_std - 1), 0.2) # standardised GARCH innovations ~ unit variance
})

test_that("dgp must be a single integer in 1:12", {
  expect_error(span_simulate(50, 2, 3, dgp = 13))
  expect_error(span_simulate(50, 2, 3, dgp = 0))
})

test_that("ncp = 0 gives near-nominal GRS size (spanning null holds)", {
  skip_on_cran()
  p <- vapply(seq_len(300), function(s) {
    set.seed(s)
    sim <- span_simulate(250, 3, 8, ncp = 0, dgp = 1)
    span_grs(sim$R1, sim$R2)$pval
  }, numeric(1))
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.03)
})

test_that("skew-t innovation path works", {
  set.seed(1)
  sim <- span_simulate(120, 2, 4, innovation = "skew-t", dynamics = "iid",
                       df = 4, xi = 0.9)
  expect_equal(dim(sim$R2), c(120L, 4L))
  expect_true(all(is.finite(sim$R2)))
})
