# =============================================================================
# Test Suite: Analytical Derivatives Validation via Numerical Differentiation
# =============================================================================
# This test file validates that all analytical gradient (gr*) and Hessian (hs*)
# functions match numerical derivatives computed via numDeriv package.
#
# Structure: 5 gradient tests + 5 Hessian tests per subfamily = 70 total tests
# Subfamilies: GKw, BKw, KKw, EKw, Mc, Kw, Beta
# =============================================================================

# -----------------------------------------------------------------------------
# Helper function to generate test parameters for each subfamily
# -----------------------------------------------------------------------------
generate_test_params <- function(subfamily, test_id) {
  # Different parameter configurations for robustness
  params_list <- list(
    gkw = list(
      list(alpha = 2.0, beta = 3.0, gamma = 1.5, delta = 0.5, lambda = 2.0),
      list(alpha = 0.8, beta = 1.2, gamma = 2.0, delta = 1.0, lambda = 0.9),
      list(alpha = 3.5, beta = 0.7, gamma = 0.8, delta = 2.0, lambda = 1.5),
      list(alpha = 1.5, beta = 2.5, gamma = 1.2, delta = 0.3, lambda = 1.8),
      list(alpha = 2.2, beta = 1.8, gamma = 0.6, delta = 1.5, lambda = 2.5)
    ),
    bkw = list(
      list(alpha = 2.0, beta = 3.0, gamma = 1.5, delta = 0.5),
      list(alpha = 0.8, beta = 1.2, gamma = 2.0, delta = 1.0),
      list(alpha = 3.5, beta = 0.7, gamma = 0.8, delta = 2.0),
      list(alpha = 1.5, beta = 2.5, gamma = 1.2, delta = 0.3),
      list(alpha = 2.2, beta = 1.8, gamma = 0.6, delta = 1.5)
    ),
    kkw = list(
      list(alpha = 2.0, beta = 3.0, delta = 0.5, lambda = 2.0),
      list(alpha = 0.8, beta = 1.2, delta = 1.0, lambda = 0.9),
      list(alpha = 3.5, beta = 0.7, delta = 2.0, lambda = 1.5),
      list(alpha = 1.5, beta = 2.5, delta = 0.3, lambda = 1.8),
      list(alpha = 2.2, beta = 1.8, delta = 1.5, lambda = 2.5)
    ),
    ekw = list(
      list(alpha = 2.0, beta = 3.0, lambda = 2.0),
      list(alpha = 0.8, beta = 1.2, lambda = 0.9),
      list(alpha = 3.5, beta = 0.7, lambda = 1.5),
      list(alpha = 1.5, beta = 2.5, lambda = 1.8),
      list(alpha = 2.2, beta = 1.8, lambda = 2.5)
    ),
    mc = list(
      list(gamma = 1.5, delta = 0.5, lambda = 2.0),
      list(gamma = 2.0, delta = 1.0, lambda = 0.9),
      list(gamma = 0.8, delta = 2.0, lambda = 1.5),
      list(gamma = 1.2, delta = 0.3, lambda = 1.8),
      list(gamma = 0.6, delta = 1.5, lambda = 2.5)
    ),
    kw = list(
      list(alpha = 2.0, beta = 3.0),
      list(alpha = 0.8, beta = 1.2),
      list(alpha = 3.5, beta = 0.7),
      list(alpha = 1.5, beta = 2.5),
      list(alpha = 2.2, beta = 1.8)
    ),
    beta = list(
      list(gamma = 2.0, delta = 3.0),
      list(gamma = 0.8, delta = 1.2),
      list(gamma = 3.5, delta = 0.7),
      list(gamma = 1.5, delta = 2.5),
      list(gamma = 2.2, delta = 1.8)
    )
  )

  params_list[[subfamily]][[test_id]]
}

# -----------------------------------------------------------------------------
# Helper function to simulate data from each subfamily
# -----------------------------------------------------------------------------
simulate_data <- function(subfamily, params, n, seed = 2203) {
  set.seed(seed)

  switch(subfamily,
    gkw = rgkw(n,
      alpha = params$alpha, beta = params$beta,
      gamma = params$gamma, delta = params$delta, lambda = params$lambda
    ),
    bkw = rbkw(n,
      alpha = params$alpha, beta = params$beta,
      gamma = params$gamma, delta = params$delta
    ),
    kkw = rkkw(n,
      alpha = params$alpha, beta = params$beta,
      delta = params$delta, lambda = params$lambda
    ),
    ekw = rekw(n, alpha = params$alpha, beta = params$beta, lambda = params$lambda),
    mc = rmc(n, gamma = params$gamma, delta = params$delta, lambda = params$lambda),
    kw = rkw(n, alpha = params$alpha, beta = params$beta),
    beta = rbeta_(n, gamma = params$gamma, delta = params$delta)
  )
}

# -----------------------------------------------------------------------------
# Helper function to get parameter vector from params list
# -----------------------------------------------------------------------------
params_to_vector <- function(subfamily, params) {
  switch(subfamily,
    gkw = c(params$alpha, params$beta, params$gamma, params$delta, params$lambda),
    bkw = c(params$alpha, params$beta, params$gamma, params$delta),
    kkw = c(params$alpha, params$beta, params$delta, params$lambda),
    ekw = c(params$alpha, params$beta, params$lambda),
    mc = c(params$gamma, params$delta, params$lambda),
    kw = c(params$alpha, params$beta),
    beta = c(params$gamma, params$delta)
  )
}

# -----------------------------------------------------------------------------
# Test tolerance (relative tolerance for numerical comparison)
# -----------------------------------------------------------------------------
GRAD_TOL <- 1e-5
HESS_TOL <- 1e-4

# =============================================================================
# GKw SUBFAMILY TESTS (5 params: alpha, beta, gamma, delta, lambda)
# =============================================================================

test_that("GKw gradient matches numerical derivative - config 1", {
  params <- generate_test_params("gkw", 1)
  data <- simulate_data("gkw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_grad <- grgkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llgkw, par_vec, data = data)


  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "GKw gradient config 1: analytical vs numerical"
  )
})

test_that("GKw gradient matches numerical derivative - config 2", {
  params <- generate_test_params("gkw", 2)
  data <- simulate_data("gkw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_grad <- grgkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llgkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "GKw gradient config 2: analytical vs numerical"
  )
})

test_that("GKw gradient matches numerical derivative - config 3", {
  params <- generate_test_params("gkw", 3)
  data <- simulate_data("gkw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_grad <- grgkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llgkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "GKw gradient config 3: analytical vs numerical"
  )
})

test_that("GKw gradient matches numerical derivative - config 4", {
  params <- generate_test_params("gkw", 4)
  data <- simulate_data("gkw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_grad <- grgkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llgkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "GKw gradient config 4: analytical vs numerical"
  )
})

test_that("GKw gradient matches numerical derivative - config 5", {
  params <- generate_test_params("gkw", 5)
  data <- simulate_data("gkw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_grad <- grgkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llgkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "GKw gradient config 5: analytical vs numerical"
  )
})

test_that("GKw Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("gkw", 1)
  data <- simulate_data("gkw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_hess <- hsgkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llgkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "GKw Hessian config 1: analytical vs numerical"
  )
})

test_that("GKw Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("gkw", 2)
  data <- simulate_data("gkw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_hess <- hsgkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llgkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "GKw Hessian config 2: analytical vs numerical"
  )
})

test_that("GKw Hessian matches numerical derivative - config 3", {
  params <- generate_test_params("gkw", 3)
  data <- simulate_data("gkw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_hess <- hsgkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llgkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "GKw Hessian config 3: analytical vs numerical"
  )
})

test_that("GKw Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("gkw", 4)
  data <- simulate_data("gkw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_hess <- hsgkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llgkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "GKw Hessian config 4: analytical vs numerical"
  )
})

test_that("GKw Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("gkw", 5)
  data <- simulate_data("gkw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("gkw", params)

  analytical_hess <- hsgkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llgkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "GKw Hessian config 5: analytical vs numerical"
  )
})

# =============================================================================
# BKw SUBFAMILY TESTS (4 params: alpha, beta, gamma, delta)
# =============================================================================

test_that("BKw gradient matches numerical derivative - config 1", {
  params <- generate_test_params("bkw", 1)
  data <- simulate_data("bkw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_grad <- grbkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llbkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "BKw gradient config 1: analytical vs numerical"
  )
})

test_that("BKw gradient matches numerical derivative - config 2", {
  params <- generate_test_params("bkw", 2)
  data <- simulate_data("bkw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_grad <- grbkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llbkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "BKw gradient config 2: analytical vs numerical"
  )
})

test_that("BKw gradient matches numerical derivative - config 3", {
  params <- generate_test_params("bkw", 3)
  data <- simulate_data("bkw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_grad <- grbkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llbkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "BKw gradient config 3: analytical vs numerical"
  )
})

test_that("BKw gradient matches numerical derivative - config 4", {
  params <- generate_test_params("bkw", 4)
  data <- simulate_data("bkw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_grad <- grbkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llbkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "BKw gradient config 4: analytical vs numerical"
  )
})

test_that("BKw gradient matches numerical derivative - config 5", {
  params <- generate_test_params("bkw", 5)
  data <- simulate_data("bkw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_grad <- grbkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llbkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "BKw gradient config 5: analytical vs numerical"
  )
})

test_that("BKw Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("bkw", 1)
  data <- simulate_data("bkw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_hess <- hsbkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "BKw Hessian config 1: analytical vs numerical"
  )
})

test_that("BKw Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("bkw", 2)
  data <- simulate_data("bkw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_hess <- hsbkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "BKw Hessian config 2: analytical vs numerical"
  )
})


test_that("BKw Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("bkw", 4)
  data <- simulate_data("bkw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_hess <- hsbkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "BKw Hessian config 4: analytical vs numerical"
  )
})

test_that("BKw Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("bkw", 5)
  data <- simulate_data("bkw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("bkw", params)

  analytical_hess <- hsbkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "BKw Hessian config 5: analytical vs numerical"
  )
})

# =============================================================================
# KKw SUBFAMILY TESTS (4 params: alpha, beta, delta, lambda)
# =============================================================================

test_that("KKw gradient matches numerical derivative - config 1", {
  params <- generate_test_params("kkw", 1)
  data <- simulate_data("kkw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_grad <- grkkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "KKw gradient config 1: analytical vs numerical"
  )
})

test_that("KKw gradient matches numerical derivative - config 2", {
  params <- generate_test_params("kkw", 2)
  data <- simulate_data("kkw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_grad <- grkkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "KKw gradient config 2: analytical vs numerical"
  )
})

test_that("KKw gradient matches numerical derivative - config 3", {
  params <- generate_test_params("kkw", 3)
  data <- simulate_data("kkw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_grad <- grkkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "KKw gradient config 3: analytical vs numerical"
  )
})

test_that("KKw gradient matches numerical derivative - config 4", {
  params <- generate_test_params("kkw", 4)
  data <- simulate_data("kkw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_grad <- grkkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "KKw gradient config 4: analytical vs numerical"
  )
})

test_that("KKw gradient matches numerical derivative - config 5", {
  params <- generate_test_params("kkw", 5)
  data <- simulate_data("kkw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_grad <- grkkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "KKw gradient config 5: analytical vs numerical"
  )
})

test_that("KKw Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("kkw", 1)
  data <- simulate_data("kkw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_hess <- hskkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "KKw Hessian config 1: analytical vs numerical"
  )
})

test_that("KKw Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("kkw", 2)
  data <- simulate_data("kkw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_hess <- hskkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "KKw Hessian config 2: analytical vs numerical"
  )
})

test_that("KKw Hessian matches numerical derivative - config 3", {
  params <- generate_test_params("kkw", 3)
  data <- simulate_data("kkw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_hess <- hskkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "KKw Hessian config 3: analytical vs numerical"
  )
})

test_that("KKw Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("kkw", 4)
  data <- simulate_data("kkw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_hess <- hskkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "KKw Hessian config 4: analytical vs numerical"
  )
})

test_that("KKw Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("kkw", 5)
  data <- simulate_data("kkw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("kkw", params)

  analytical_hess <- hskkw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "KKw Hessian config 5: analytical vs numerical"
  )
})

# =============================================================================
# EKw SUBFAMILY TESTS (3 params: alpha, beta, lambda)
# =============================================================================

test_that("EKw gradient matches numerical derivative - config 1", {
  params <- generate_test_params("ekw", 1)
  data <- simulate_data("ekw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_grad <- grekw(par_vec, data)
  numerical_grad <- numDeriv::grad(llekw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "EKw gradient config 1: analytical vs numerical"
  )
})

test_that("EKw gradient matches numerical derivative - config 2", {
  params <- generate_test_params("ekw", 2)
  data <- simulate_data("ekw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_grad <- grekw(par_vec, data)
  numerical_grad <- numDeriv::grad(llekw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "EKw gradient config 2: analytical vs numerical"
  )
})

test_that("EKw gradient matches numerical derivative - config 3", {
  params <- generate_test_params("ekw", 3)
  data <- simulate_data("ekw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_grad <- grekw(par_vec, data)
  numerical_grad <- numDeriv::grad(llekw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "EKw gradient config 3: analytical vs numerical"
  )
})

test_that("EKw gradient matches numerical derivative - config 4", {
  params <- generate_test_params("ekw", 4)
  data <- simulate_data("ekw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_grad <- grekw(par_vec, data)
  numerical_grad <- numDeriv::grad(llekw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "EKw gradient config 4: analytical vs numerical"
  )
})

test_that("EKw gradient matches numerical derivative - config 5", {
  params <- generate_test_params("ekw", 5)
  data <- simulate_data("ekw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_grad <- grekw(par_vec, data)
  numerical_grad <- numDeriv::grad(llekw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "EKw gradient config 5: analytical vs numerical"
  )
})

test_that("EKw Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("ekw", 1)
  data <- simulate_data("ekw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_hess <- hsekw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llekw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "EKw Hessian config 1: analytical vs numerical"
  )
})

test_that("EKw Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("ekw", 2)
  data <- simulate_data("ekw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_hess <- hsekw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llekw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "EKw Hessian config 2: analytical vs numerical"
  )
})

test_that("EKw Hessian matches numerical derivative - config 3", {
  params <- generate_test_params("ekw", 3)
  data <- simulate_data("ekw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_hess <- hsekw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llekw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "EKw Hessian config 3: analytical vs numerical"
  )
})

test_that("EKw Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("ekw", 4)
  data <- simulate_data("ekw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_hess <- hsekw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llekw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "EKw Hessian config 4: analytical vs numerical"
  )
})

test_that("EKw Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("ekw", 5)
  data <- simulate_data("ekw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("ekw", params)

  analytical_hess <- hsekw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llekw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "EKw Hessian config 5: analytical vs numerical"
  )
})

# =============================================================================
# Mc SUBFAMILY TESTS (3 params: gamma, delta, lambda)
# =============================================================================

test_that("Mc gradient matches numerical derivative - config 1", {
  params <- generate_test_params("mc", 1)
  data <- simulate_data("mc", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_grad <- grmc(par_vec, data)
  numerical_grad <- numDeriv::grad(llmc, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Mc gradient config 1: analytical vs numerical"
  )
})

test_that("Mc gradient matches numerical derivative - config 2", {
  params <- generate_test_params("mc", 2)
  data <- simulate_data("mc", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_grad <- grmc(par_vec, data)
  numerical_grad <- numDeriv::grad(llmc, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Mc gradient config 2: analytical vs numerical"
  )
})

test_that("Mc gradient matches numerical derivative - config 3", {
  params <- generate_test_params("mc", 3)
  data <- simulate_data("mc", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_grad <- grmc(par_vec, data)
  numerical_grad <- numDeriv::grad(llmc, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Mc gradient config 3: analytical vs numerical"
  )
})

test_that("Mc gradient matches numerical derivative - config 4", {
  params <- generate_test_params("mc", 4)
  data <- simulate_data("mc", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_grad <- grmc(par_vec, data)
  numerical_grad <- numDeriv::grad(llmc, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Mc gradient config 4: analytical vs numerical"
  )
})

test_that("Mc gradient matches numerical derivative - config 5", {
  params <- generate_test_params("mc", 5)
  data <- simulate_data("mc", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_grad <- grmc(par_vec, data)
  numerical_grad <- numDeriv::grad(llmc, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Mc gradient config 5: analytical vs numerical"
  )
})

test_that("Mc Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("mc", 1)
  data <- simulate_data("mc", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_hess <- hsmc(par_vec, data)
  numerical_hess <- numDeriv::hessian(llmc, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Mc Hessian config 1: analytical vs numerical"
  )
})

test_that("Mc Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("mc", 2)
  data <- simulate_data("mc", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_hess <- hsmc(par_vec, data)
  numerical_hess <- numDeriv::hessian(llmc, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Mc Hessian config 2: analytical vs numerical"
  )
})

test_that("Mc Hessian matches numerical derivative - config 3", {
  params <- generate_test_params("mc", 3)
  data <- simulate_data("mc", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_hess <- hsmc(par_vec, data)
  numerical_hess <- numDeriv::hessian(llmc, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Mc Hessian config 3: analytical vs numerical"
  )
})

test_that("Mc Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("mc", 4)
  data <- simulate_data("mc", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_hess <- hsmc(par_vec, data)
  numerical_hess <- numDeriv::hessian(llmc, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Mc Hessian config 4: analytical vs numerical"
  )
})

test_that("Mc Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("mc", 5)
  data <- simulate_data("mc", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("mc", params)

  analytical_hess <- hsmc(par_vec, data)
  numerical_hess <- numDeriv::hessian(llmc, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Mc Hessian config 5: analytical vs numerical"
  )
})

# =============================================================================
# Kw SUBFAMILY TESTS (2 params: alpha, beta)
# =============================================================================

test_that("Kw gradient matches numerical derivative - config 1", {
  params <- generate_test_params("kw", 1)
  data <- simulate_data("kw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_grad <- grkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Kw gradient config 1: analytical vs numerical"
  )
})

test_that("Kw gradient matches numerical derivative - config 2", {
  params <- generate_test_params("kw", 2)
  data <- simulate_data("kw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_grad <- grkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Kw gradient config 2: analytical vs numerical"
  )
})

test_that("Kw gradient matches numerical derivative - config 3", {
  params <- generate_test_params("kw", 3)
  data <- simulate_data("kw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_grad <- grkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Kw gradient config 3: analytical vs numerical"
  )
})

test_that("Kw gradient matches numerical derivative - config 4", {
  params <- generate_test_params("kw", 4)
  data <- simulate_data("kw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_grad <- grkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Kw gradient config 4: analytical vs numerical"
  )
})

test_that("Kw gradient matches numerical derivative - config 5", {
  params <- generate_test_params("kw", 5)
  data <- simulate_data("kw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_grad <- grkw(par_vec, data)
  numerical_grad <- numDeriv::grad(llkw, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Kw gradient config 5: analytical vs numerical"
  )
})

test_that("Kw Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("kw", 1)
  data <- simulate_data("kw", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_hess <- hskw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Kw Hessian config 1: analytical vs numerical"
  )
})

test_that("Kw Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("kw", 2)
  data <- simulate_data("kw", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_hess <- hskw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Kw Hessian config 2: analytical vs numerical"
  )
})

test_that("Kw Hessian matches numerical derivative - config 3", {
  params <- generate_test_params("kw", 3)
  data <- simulate_data("kw", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_hess <- hskw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Kw Hessian config 3: analytical vs numerical"
  )
})

test_that("Kw Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("kw", 4)
  data <- simulate_data("kw", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_hess <- hskw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Kw Hessian config 4: analytical vs numerical"
  )
})

test_that("Kw Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("kw", 5)
  data <- simulate_data("kw", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("kw", params)

  analytical_hess <- hskw(par_vec, data)
  numerical_hess <- numDeriv::hessian(llkw, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Kw Hessian config 5: analytical vs numerical"
  )
})

# =============================================================================
# Beta SUBFAMILY TESTS (2 params: gamma, delta)
# =============================================================================

test_that("Beta gradient matches numerical derivative - config 1", {
  params <- generate_test_params("beta", 1)
  data <- simulate_data("beta", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_grad <- grbeta(par_vec, data)
  numerical_grad <- numDeriv::grad(llbeta, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Beta gradient config 1: analytical vs numerical"
  )
})

test_that("Beta gradient matches numerical derivative - config 2", {
  params <- generate_test_params("beta", 2)
  data <- simulate_data("beta", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_grad <- grbeta(par_vec, data)
  numerical_grad <- numDeriv::grad(llbeta, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Beta gradient config 2: analytical vs numerical"
  )
})

test_that("Beta gradient matches numerical derivative - config 3", {
  params <- generate_test_params("beta", 3)
  data <- simulate_data("beta", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_grad <- grbeta(par_vec, data)
  numerical_grad <- numDeriv::grad(llbeta, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Beta gradient config 3: analytical vs numerical"
  )
})

test_that("Beta gradient matches numerical derivative - config 4", {
  params <- generate_test_params("beta", 4)
  data <- simulate_data("beta", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_grad <- grbeta(par_vec, data)
  numerical_grad <- numDeriv::grad(llbeta, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Beta gradient config 4: analytical vs numerical"
  )
})

test_that("Beta gradient matches numerical derivative - config 5", {
  params <- generate_test_params("beta", 5)
  data <- simulate_data("beta", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_grad <- grbeta(par_vec, data)
  numerical_grad <- numDeriv::grad(llbeta, par_vec, data = data)

  expect_equal(analytical_grad, numerical_grad,
    tolerance = GRAD_TOL,
    info = "Beta gradient config 5: analytical vs numerical"
  )
})

test_that("Beta Hessian matches numerical derivative - config 1", {
  params <- generate_test_params("beta", 1)
  data <- simulate_data("beta", params, n = 1500, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_hess <- hsbeta(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbeta, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Beta Hessian config 1: analytical vs numerical"
  )
})

test_that("Beta Hessian matches numerical derivative - config 2", {
  params <- generate_test_params("beta", 2)
  data <- simulate_data("beta", params, n = 1200, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_hess <- hsbeta(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbeta, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Beta Hessian config 2: analytical vs numerical"
  )
})

test_that("Beta Hessian matches numerical derivative - config 3", {
  params <- generate_test_params("beta", 3)
  data <- simulate_data("beta", params, n = 1800, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_hess <- hsbeta(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbeta, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Beta Hessian config 3: analytical vs numerical"
  )
})

test_that("Beta Hessian matches numerical derivative - config 4", {
  params <- generate_test_params("beta", 4)
  data <- simulate_data("beta", params, n = 1000, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_hess <- hsbeta(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbeta, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Beta Hessian config 4: analytical vs numerical"
  )
})

test_that("Beta Hessian matches numerical derivative - config 5", {
  params <- generate_test_params("beta", 5)
  data <- simulate_data("beta", params, n = 2000, seed = 2203)
  par_vec <- params_to_vector("beta", params)

  analytical_hess <- hsbeta(par_vec, data)
  numerical_hess <- numDeriv::hessian(llbeta, par_vec, data = data)

  expect_equal(analytical_hess, numerical_hess,
    tolerance = HESS_TOL,
    info = "Beta Hessian config 5: analytical vs numerical"
  )
})
