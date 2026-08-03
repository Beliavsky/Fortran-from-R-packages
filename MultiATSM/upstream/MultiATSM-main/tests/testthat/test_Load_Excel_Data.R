library(testthat)
library(MultiATSM)


# Unit test for Load_Excel_Data
test_that("Load_Excel_Data loads package example files correctly", {
  # Test with actual package data files
  macro_file <- system.file("extdata", "MacroData.xlsx", package = "MultiATSM")
  yields_file <- system.file("extdata", "YieldsData.xlsx", package = "MultiATSM")

  # Test that files exist (sanity check)
  expect_true(file.exists(macro_file))
  expect_true(file.exists(yields_file))

  # Test function with package files
  expect_silent(macro_data <- Load_Excel_Data(macro_file))
  expect_silent(yields_data <- Load_Excel_Data(yields_file))

  # Test basic structure
  expect_type(macro_data, "list")
  expect_type(yields_data, "list")
  expect_true(length(macro_data) > 0)
  expect_true(length(yields_data) > 0)
  expect_named(macro_data)
  expect_named(yields_data)
})

test_that("Load_Excel_Data handles invalid inputs", {
  # Test with non-existent file
  expect_error(Load_Excel_Data("nonexistent_file.xlsx"))

  # Test with NULL input
  expect_error(Load_Excel_Data(NULL))

  # Test with non-Excel file
  temp_file <- tempfile(fileext = ".txt")
  writeLines("test", temp_file)
  expect_error(Load_Excel_Data(temp_file))
  unlink(temp_file)
})

test_that("Load_Excel_Data returns correct structure", {
  macro_file <- system.file("extdata", "MacroData.xlsx", package = "MultiATSM")
  skip_if_not(file.exists(macro_file))

  result <- Load_Excel_Data(macro_file)

  # Test generic structure (without specific content checks)
  expect_type(result, "list")
  expect_true(length(result) > 0)
  expect_true(all(sapply(result, is.data.frame)))
  expect_named(result)
})

test_that("Load_Excel_Data preserves sheet names and order", {
  # Use an existing file from your package that has multiple sheets
  test_file <- system.file("extdata", "MacroData.xlsx", package = "MultiATSM")
  skip_if_not(file.exists(test_file))

  # First, check what sheets exist in the actual file
  actual_sheets <- readxl::excel_sheets(test_file)
  skip_if(length(actual_sheets) < 2) # Skip if file doesn't have multiple sheets

  result <- Load_Excel_Data(test_file)

  # Test that sheet names and order are preserved from the actual file
  expect_equal(names(result), actual_sheets)
  expect_true(all(sapply(result, is.data.frame)))
})


test_that("Load_Excel_Data integrates with package workflow", {
  # Test that the function works with package's expected data flow
  macro_file <- system.file("extdata", "MacroData.xlsx", package = "MultiATSM")

  # Load data using function
  macro_data <- Load_Excel_Data(macro_file)

  # Test that data is suitable for downstream package functions
  expect_true(all(sapply(macro_data, is.data.frame)))
  expect_true(all(sapply(macro_data, function(df) nrow(df) > 0)))

  # Test that data has expected financial/time series structure
  first_sheet <- macro_data[[1]]
  expect_true(any(grepl("date|Date|TIME|Time", names(first_sheet), ignore.case = TRUE)) ||
    nrow(first_sheet) > 1) # Either has date column or is time series data
})
