test_that("printing percent", {

  format_percent(0.52366) |> expect_snapshot()
  format_percent(0.52366, accuracy = 0.01) |> expect_snapshot()
  format_percent(list(a = 0.52366, b = 0.23456, c = "test")) |> 
    expect_snapshot()
  format_percent(
    list(
      a = 0.52366, 
      b = 0.23456, 
      c = "test", 
      d = list(
        a = 0.52366, 
        b = 0.23456, 
        c = "test"
      )), 
    accuracy = 0.01
  ) |> 
    expect_snapshot()
})

test_that("printing currency", {

  format_currency(234) |> expect_snapshot()
  format_currency(234, prefix = "$") |> expect_snapshot()
  format_currency(234, suffix = " PLN") |> expect_snapshot()
  format_currency(1234567.123456) |> expect_snapshot()
  format_currency(1234567.123456, accuracy = 0.01) |> expect_snapshot()
  format_currency(list(a = 234, b = 1234567.123456, c = "test")) |> 
    expect_snapshot()
  format_currency(
    list(
      a = 234, 
      b = 1234567.123456, 
      c = "test", 
      d = list(
        a = 234, 
        b = 1234567.123456, 
        c = "test"
      )), 
    accuracy = 0.01
  ) |> 
    expect_snapshot()
})

test_that("printing currency with padding", {

  expect_equal(
    format_currency(234),
    "234"
  )
  expect_equal(
    format_currency(234, min_length = 6),
    "   234"
  )
  expect_equal(
    format_currency(2343, min_length = 6),
    " 2,343"
  )
  expect_equal(
    format_currency(2343, min_length = 7),
    "  2,343"
  )
  expect_equal(
    format_currency(0, min_length = 6),
    "     0"
  )
})

test_that("getting default current date", {

  withr::local_options(
    R4GPF.current_date = NULL
  )

  expect_equal(
    get_current_date(),
    Sys.Date()
  )
})

test_that("getting custom current date", {

  test_current_date <- "2020-01-01"

  withr::local_options(
    R4GPF.current_date = test_current_date
  )

  expect_equal(
    get_current_date(),
    lubridate::as_date(test_current_date)
  )
})

test_that("normalizing value", {
  
  expect_snapshot(normalize(1:10))
  expect_snapshot(normalize(1:10, 10, 20))
})

test_that("generating multiple random seeds", {

  withr::local_seed(1234)

  seeds <- generate_random_seed_vector(5) 
  if (interactive()) print(seeds)
  
  expect_true(is.integer(seeds) |> all())

  expect_equal(
    length(seeds), 
    5
  )

  expect_equal(
    unique(seeds) |> length(),
    5
  )
})
