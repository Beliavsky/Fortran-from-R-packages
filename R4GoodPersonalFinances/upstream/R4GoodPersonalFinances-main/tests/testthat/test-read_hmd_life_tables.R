test_that("reading hmd life tables for all genders", {

  temp_dir <- withr::local_tempdir()

  path <- testthat::test_path("fixtures", "life-tables", "usa")
  
  fs::dir_copy(
    path      = path,
    new_path  = temp_dir,
    overwrite = TRUE
  )

  test_life_tables <- read_hmd_life_tables(path = temp_dir) 

  expect_true(is.character(test_life_tables$sex))
  expect_true(is.integer(test_life_tables$age))
  expect_true(is.integer(test_life_tables$year))
  expect_true(is.numeric(test_life_tables$mortality_rate))
  expect_true(is.numeric(test_life_tables$life_expectancy))
})

test_that("reading hmd life tables for selected genders", {

  temp_dir <- withr::local_tempdir()

  path <- testthat::test_path("fixtures", "life-tables", "usa")

  fs::dir_copy(
    path      = path,
    new_path  = temp_dir,
    overwrite = TRUE
  )

  fs::file_delete(fs::path(temp_dir, "bltper_1x1.txt"))

  expect_warning(
    test_life_tables <- read_hmd_life_tables(path = temp_dir) 
  )

  expect_true(is.data.frame(test_life_tables))
})
