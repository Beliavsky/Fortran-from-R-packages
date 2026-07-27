test_that("str_length of factor is length of level", {
  expect_equal(.diag_of(matrix(c(1,2,3,3,2,1,2,3,1), ncol = 3, nrow = 3)), c(1,2,1))
})
