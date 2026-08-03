test_that("prepare_data filters and processes data correctly", {
  # Ensure the required packages are loaded for the test
  library(lubridate)

  # Example input data
  data <- data.frame(
    Date = c("19990101", "19990115", "19990201", "19990301", "19990315", "19990401"),
    Var1 = c(1, 2, -99.99, 4, 5, -99.99),
    Var2 = c(3, -99.99, 6, 7, 8, 9),
    Var3 = c(10, 11, 12, 13, -99.99, 15)
  )

  # Run the function with default date range
  result <- prepare_data(data, date_column_index = 1, start_date = '19990101', end_date = '19990430')

  # Check the structure of the result
  expect_type(result, "list")
  expect_true("x" %in% names(result))
  expect_true("mon" %in% names(result))
  expect_true("count" %in% names(result))
  expect_true("Date" %in% names(result))
})
