test_that("span_gl_a returns correct structure and values", {
  set.seed(1234)
  R1 <- matrix(rnorm(300), 100, 3)
  R2 <- matrix(rnorm(200), 100, 2)

  res <- span_gl_a(R1, R2, control = list(totsim = 100, do_trace = FALSE))

  expect_type(res, "list")
  expect_named(res, c("pval_LMC", "pval_BMC", "stat", "Decisions", "Decisions_string", "H0"))
  expect_true(is.numeric(res$pval_LMC) && res$pval_LMC >= 0 && res$pval_LMC <= 1)
  expect_true(is.numeric(res$pval_BMC) && res$pval_BMC >= 0 && res$pval_BMC <= 1)
  expect_true(res$Decisions_string %in% c("Accept", "Reject", "Inconclusive"))
  expect_equal(res$H0, "alpha = 0")
})
