library(testthat)
library(MultiATSM)


test_that("LoadData loads correct datasets for each paper", {
  # Test all valid cases
  valid_cases <- c("BR_2017", "CM_2024", "CM_2023")

  for (case in valid_cases) {
    expect_silent(LoadData(case))
  }
})
