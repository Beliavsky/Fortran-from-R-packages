test_that("model construction rejects invalid parameters", {
  expect_error(varmapack_model(Sig = 1), "square matrix")
  expect_error(varmapack_model(Sig = matrix(c(1, 2, 0, 1), 2, 2)),
               "symmetric")
  expect_error(varmapack_model(A = matrix(1, 2, 2), Sig = diag(1)),
               "dimensions")
  expect_error(varmapack_model(B = array(NA_real_, c(1, 1, 1)), Sig = diag(1)),
               "finite")
  expect_error(varmapack_model(C = matrix(1, 2, 1), Sig = diag(1)),
               "dimensions")
  expect_error(varmapack_model(Sig = diag(2), mu = 1:3), "mu")
  expect_error(varmapack_model(C = matrix(1, 1, 1), Sig = diag(1), mu = 0),
               "not supported")
})

test_that("VARMA simulation rejects invalid arguments", {
  model <- varmapack_model(A = matrix(0.4, 1, 1), Sig = matrix(1, 1, 1))
  expect_error(model$sim(0), "length")
  expect_error(model$sim(2.5), "length")
  expect_error(model$sim(4, nrep = 0), "nrep")
  expect_error(model$sim(4, return_shocks = 1), "return_shocks")
  expect_error(model$sim(4, return_shocks = "x"), "return_shocks")
  expect_error(model$sim(4, return_shocks = c(TRUE, FALSE)), "return_shocks")
  expect_error(model$sim(4, rng = new.env()), "rng")
  expect_error(model$sim(4, X0 = 1:2), "X0")
  expect_error(model$sim(4, nrep = 2, X0 = array(0, c(1, 1, 3))), "X0")
  expect_error(model$sim(4, z = matrix(0, 1, 4)), "VARMAX")
  unstable <- varmapack_model(A = matrix(1.1, 1, 1), Sig = matrix(1, 1, 1))
  expect_error(unstable$sim(4), "nonstationary")
})

test_that("analysis functions reject invalid arguments", {
  model <- varmapack_model(Sig = diag(2))
  expect_error(model$acvf(-1), "maxlag")
  expect_error(model$psi(-1), "maxlag")
  expect_error(model$irf(-1), "maxlag")
  expect_error(varmapack_autocov(1:4, 1), "matrix")
  expect_error(varmapack_autocov(matrix(c(1, NA), 1, 2), 1), "finite")
  expect_error(varmapack_autocov(matrix(1:4, 1, 4), 4), "maxlag")
  expect_error(varmapack_autocov(matrix(1:4, 1, 4), 1, "bad"), "norm")
  expect_error(varmapack_cov2corr(1:4), "matrix")
  expect_error(varmapack_cov2corr(matrix(1, 2, 3)), "dimensions")
  expect_error(varmapack_cov2corr(matrix(c(0, 0, 0, 1), 2, 2)),
               "nonpositive variance")
})

test_that("VARMAX simulation rejects missing and mismatched paths", {
  model <- varmapack_model(A = matrix(0.2, 1, 1), C = matrix(0.3, 1, 1),
                           Sig = matrix(1, 1, 1))
  expect_error(model$sim(4, X0 = matrix(0, 1, 1)), "z")
  expect_error(model$sim(4, z = matrix(0, 1, 4)), "X0")
  expect_error(model$sim(4, X0 = matrix(0, 1, 1), z = matrix(0, 2, 4)), "z")
  expect_error(model$sim(4, nrep = 2, X0 = array(0, c(1, 1, 3)),
                         z = matrix(0, 1, 4)), "X0")
  expect_error(model$sim(4, nrep = 2, X0 = matrix(0, 1, 1),
                         z = array(0, c(1, 4, 3))), "z")
})

test_that("testcase construction rejects invalid arguments", {
  expect_error(varmapack_testcase("does-not-exist"), "unknown testcase")
  expect_error(varmapack_testcase(list()), "which")
  expect_error(varmapack_testcase("rho"), "p")
  expect_error(varmapack_testcase(.Machine$integer.max + 1), "which")
  expect_error(varmapack_testcase("rho", p = 1, q = 1, r = 1, rho = -1),
               "rho")
  expect_error(varmapack_testcase("random", p = 1, q = 1, r = 1,
                                  rng = new.env()), "rng")
})
