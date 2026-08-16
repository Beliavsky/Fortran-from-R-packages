dslnex <- function(x) {
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1]-1) + x[2]^3 - 2
  y
}

jacdsln <- function(x) {
  n <- length(x)
  Df <- matrix(numeric(n*n),n,n)
  Df[1,1] <- 2*x[1]
  Df[1,2] <- 2*x[2]
  Df[2,1] <- exp(x[1]-1)
  Df[2,2] <- 3*x[2]^2

  Df
}

test_that("nleqslv error conditions work", {
  xstart <- c(1, 2, 3)
  z <- nleqslv(xstart, common_test_f, common_test_jac, method="Newton",
               control=list(trace=0, allowSingular=TRUE))
  expect_equal(z$x, common_test_xsol, tolerance = 1E-6)

  expect_error(z <- nleqslv(xstart, common_test_f, "Jac", method="Newton"))

  expect_error(z <- nleqslv(xstart, common_test_f, common_test_jac, method="Newton",
                            control=list(delta = as.Date("2026-01-01"))))

  xstart <- c(2,0.5)
  expect_output(
    expect_error(z <- nleqslv(xstart, dslnex, jacobian = TRUE, method="Newton",
                            global="none", xscalm="fixed",
                            control=list(trace=1))))
  expect_output(
    expect_error(z <- nleqslv(xstart, dslnex, jacobian = TRUE, method="Newton",
                            global="none", xscalm="fixed",
                            control=list(trace=1))))

  expect_error(nleqslv(xstart, dslnex, control=list(xtol=Inf)))
  expect_error(nleqslv(xstart, dslnex, control=list(ftol=Inf)))
  expect_error(nleqslv(xstart, dslnex, control=list(btol=Inf)))
  expect_error(nleqslv(xstart, dslnex, control=list(sigma=Inf)))
  expect_error(nleqslv(xstart, dslnex, control=list(stepmax=Inf)))
  expect_error(nleqslv(xstart, dslnex, control=list(delta=Inf)))
  expect_error(nleqslv(xstart, dslnex, control=list(cndtol=Inf)))

  expect_error(nleqslv(LETTERS[1:2], dslnex))

  expect_error(nleqslv(c(0, Inf), dslnex))

  expect_error(nleqslv(c(1.0, 2.0), function(x) return(LETTERS[1:2])))

  expect_error(nleqslv(c(1.0, 2.0), function(x) return(c(0, Inf))))

  expect_error(nleqslv(c(1.0, 2.0), function(x) numeric(3)))

  expect_error(nleqslv(c(1.0, 2.0), dslnex, control=list(scalex=c(1,2,3))))

  expect_error(nleqslv(c(1.0, 2.0), dslnex, control=list(scalex=c(1, Inf))))
})

test_that("nleqslv - test all combinations of parameters", {
  xstart <- c(2,0.5)

  for (m in c("Newton", "Broyden")) {
    for (g in c("dbldog", "pwldog",
                "cline", "qline", "gline", "hook", "none")) {
      for (x in c("fixed", "auto")) {
        if (!(m == "Newton" & g == "none")) {
          expect_output(nleqslv(xstart, dslnex, jacobian = TRUE, method=m,
                                  global=g, xscalm=x,
                                  control=list(trace=1)))
        }
      }
    }
  }
})
