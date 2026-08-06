test_that("f_getpv returns named vector of p-values", {
  set.seed(1)
  x <- matrix(rnorm(300), nrow = 100, ncol = 3)
  u <- rnorm(100)

  pvals <- spantest:::f_getpv(u, x, ks = c(1/3), L = c(0, 2))

  expect_type(pvals, "double")
  expect_true(all(pvals >= 0 & pvals <= 1))
  expect_true(all(grepl("^CCT[ad]{1,2}_L\\d+_k\\d+$", names(pvals))))
})
