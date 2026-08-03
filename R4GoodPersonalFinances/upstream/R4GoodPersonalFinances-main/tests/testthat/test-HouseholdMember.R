test_that("setting birth date", {

  test_birth_date <- "1980-07-15"
  
  member <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  
  expect_true(
    inherits(member, "HouseholdMember")
  )
  
  expect_true(
    inherits(member$get_birth_date(), "Date")
  )
  
  expect_equal(
    as.character(member$get_birth_date()), 
    test_birth_date
  )

  expect_equal(
    member$get_name(),
    "test_name"
  )
})

test_that("setting max age", {
  
  test_birth_date <- "1980-07-15"
  test_max_age    <- 120

  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  expect_equal(members$max_age, 100)

  members$max_age <- test_max_age
  expect_equal(
    members$max_age,
    test_max_age
  )
})

test_that("calculating age", {

  test_birth_date   <- "1980-07-15"
  test_current_date <- 
    c("2020-07-15", "2025-01-01", "2080-07-15", 
      "2080-07-16", "2081-07-11", "2081-07-18")

  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  expect_equal(members$max_age, 100)
  
  expect_snapshot_value(
    style = "json2",
    members$calc_age(current_date = test_current_date)
  )
})

test_that("calculating max lifespan", {
  
  test_birth_date   <- "1980-07-15"
  test_current_date <- 
    c("2020-07-15", "2025-01-01", "2080-07-15", "2080-07-16")
  
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )

  expect_snapshot_value(
    style = "json2",
    members$get_lifespan(current_date = test_current_date)
  )
})

test_that("setting gompertz parameters", {
  
  test_birth_date   <- "1980-07-15"
  test_current_date <- "2020-07-15"
  
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  expect_null(members$mode)
  expect_null(members$dispersion)

  members$mode <- 88
  expect_equal(members$mode, 88)

  members$dispersion <- 10
  expect_equal(members$dispersion, 10)
})

test_that("calculating gompertz survival probability", {

  test_birth_date   <- "1955-07-15"
  test_current_date <- "2020-07-15"
  
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  members$mode       <- 80
  members$dispersion <- 10
  
  expect_equal(
    members$calc_survival_probability(
      target_age   = 85,
      current_date = test_current_date
    ),
    0.2404, 
    tolerance = 0.001
  )
})

test_that("setting age event", {
  
  test_birth_date <- "1955-07-15"
  
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )

  expect_true(is.list(members$get_events()))
  expect_equal(NROW(members$get_events()), 0)
  
  members$set_event("retirement", 65)
  expect_equal(members$get_events()$retirement$start_age, 65)
  expect_equal(members$get_events()$retirement$end_age, Inf)
  
  members$set_event("social_security", 70)
  expect_equal(members$get_events()$social_security$start_age, 70)
  expect_equal(members$get_events()$social_security$end_age, Inf)
  
  members$set_event("kid", 20, years = 20)
  expect_equal(members$get_events()$kid$start_age, 20)
  expect_equal(members$get_events()$kid$end_age, 40 - 1)
  
  expect_equal(NROW(members$get_events()), 3)
})

test_that("cloning works", {

  test_birth_date <- "1955-07-15"
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  
  members$set_event("retirement", 65)
  expect_equal(members$get_events()$retirement$start_age, 65)
  
  cloned_hm <- members$clone(deep = TRUE)
  cloned_hm$set_event("retirement", 100)
  expect_equal(cloned_hm$get_events()$retirement$start_age, 100)
})

test_that("calculating life expectancy", {

  member <- HouseholdMember$new(
    name       = "Isabela",
    birth_date = Sys.Date() - lubridate::years(25),
    mode       = 91,
    dispersion = 8.88
  )

  expect_equal(
    round(member$calc_life_expectancy()), 
    86
  )
})
