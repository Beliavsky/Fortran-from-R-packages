test_that("test HW adagio", {
  if (requireNamespace("adagio", quietly = TRUE)) {

    n <- 10
    fn <- adagio::fnRosenbrock
    gr <- adagio::grRosenbrock
    lb <- c(rep(-2.5, n)); ub = -lb
    x0 <- rep(0.5, n)

    sol1 <- lbfgsb3c(x0, fn, gr = gr, lower = lb, upper = ub)


    expect_true(sol1$value < 1e-4)

  }
})
