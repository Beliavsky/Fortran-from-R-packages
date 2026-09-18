test_that("named and indexed testcases create models", {
  named <- varmapack_testcase("smallARMA1")
  indexed <- varmapack_testcase(8)
  expect_equal(named$A, indexed$A)
  expect_equal(named$B, indexed$B)
  expect_equal(named$Sig, indexed$Sig)
})

test_that("unnamed testcases support rho and Randompack RNGs", {
  rho_case <- varmapack_testcase("rho", p = 2, q = 1, r = 2, rho = 0.8)
  expect_equal(rho_case$specrad(), 0.8, tolerance = 2e-6)
  rng1 <- randompack::randompack_rng()
  rng2 <- randompack::randompack_rng()
  rng1$seed(321)
  rng2$seed(321)
  random1 <- varmapack_testcase("random", p = 2, q = 1, r = 2, rng = rng1)
  random2 <- varmapack_testcase("random", p = 2, q = 1, r = 2, rng = rng2)
  expect_identical(random1$A, random2$A)
  expect_identical(random1$B, random2$B)
  expect_identical(random1$Sig, random2$Sig)
  deterministic <- varmapack_testcase("deterministic", p = 2, q = 1, r = 2)
  expect_equal(dim(deterministic$A), c(2L, 2L, 2L))
  expect_equal(dim(deterministic$B), c(2L, 2L, 1L))
})

test_that("testcase listing has all named cases", {
  cases <- varmapack_testcases()
  expect_equal(nrow(cases), 16L)
  expect_named(cases, c("index", "name", "p", "q", "r"))
  tiny <- cases[cases$name == "tinyAR", ]
  large_ar <- cases[cases$name == "largeAR", ]
  large_arma <- cases[cases$name == "largeARMA", ]
  expect_equal(unname(unlist(tiny[c("index", "p", "q", "r")])),
               c(1L, 1L, 0L, 1L))
  expect_equal(unname(unlist(large_ar[c("index", "p", "q", "r")])),
               c(15L, 5L, 0L, 7L))
  expect_equal(unname(unlist(large_arma[c("index", "p", "q", "r")])),
               c(16L, 3L, 3L, 7L))
})
