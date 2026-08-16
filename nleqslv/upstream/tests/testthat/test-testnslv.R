test_that("testnslv works", {
  xstart <- c(1, 2, 3)
  expect_no_error(testnslv(xstart, common_test_f, common_test_jac,
                           method="Newton", Nrep = 0,
                            control = list(allowSingular=TRUE)))
  expect_no_error(testnslv(xstart, common_test_f, common_test_jac,
                           method="Newton", Nrep = 1,
                           control = list(allowSingular=TRUE)))
  expect_no_error(testnslv(xstart, common_test_f, common_test_jac,
                           method="Newton", Nrep = 10,
                           control = list(allowSingular=TRUE)))

  expect_error(z <- nleqslv(xstart, common_test_f, "Jac", method="Newton"))

  expect_error(z <- nleqslv(xstart, common_test_f, common_test_jac, method="Newton",
                            control=list(delta = as.Date("2026-01-01"))))

})

test_that("print.testnslv", {
  xstart <- c(1, 2, 3)
  z <- testnslv(xstart, common_test_f, common_test_jac,
                           method="Newton", Nrep = 0,
                           control = list(allowSingular=TRUE))
  expect_output(print(z))

  expect_error(print.nleqlsv("A"))

  expect_no_error(testnslv(xstart, common_test_f, common_test_jac,
                          method="Newton", Nrep = 0,
                          control = list(allowSingular=TRUE)))
  z$title <- "Test Title"
  expect_output(print(z))
  z$title <- NULL
  expect_output(print(z))
})
