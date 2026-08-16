test_that("seachZeros works", {
  xstart <- matrix(c(1,2,3,
                     4,5,6), nrow = 2, ncol = 3, byrow = TRUE)
  z <- searchZeros(xstart, common_test_f)

  expect_error(searchZeros("A", common_test_f))
  expect_error(searchZeros(1:3, common_test_f))
  expect_no_error(searchZeros(xstart, common_test_f, digits = as.numeric(NA)))
  expect_error(searchZeros(xstart, common_test_f, digits = "A"))
  colnames(xstart) <- LETTERS[1:3]
  expect_no_error(searchZeros(xstart, common_test_f))

  xstart[2,2] <- NA
  expect_error(searchZeros(xstart, common_test_f))

  expect_error(searchZeros(1:3, common_test_f))

  expect_error(searchZeros(matrix(nrow = 0, ncol = 5), common_test_f))
})
