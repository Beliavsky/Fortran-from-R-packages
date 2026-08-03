test_that("adding household members", {
  
  h <- Household$new()

  expect_true(is.list(h$get_members()))
  expect_true(length(h$get_members()) == 0)
  
  h$add_member(
    HouseholdMember$new(
      name       = "test_member_1",  
      birth_date = "1980-07-15",
        mode       = 80,
        dispersion = 10
    )
  )
  expect_equal(
    length(h$get_members()), 
    1
  )
  expect_true(
    inherits(h$get_members()[["test_member_1"]], "HouseholdMember")
  )
  
  expect_error(
    h$add_member(
      HouseholdMember$new(
        name       = "test_member_1",  
        birth_date = "1980-07-15"
      )
    )
  )

  h$add_member(
    HouseholdMember$new(
      name       = "test_member_2",  
      birth_date = "1980-07-15"
    )
  )
  expect_equal(
    length(h$get_members()), 
    2
  )
  expect_true(
    inherits(h$get_members()[["test_member_1"]], "HouseholdMember")
  )
  expect_true(
    inherits(h$get_members()[["test_member_2"]], "HouseholdMember")
  )

  if (interactive()) print(h)
})

test_that("calculating household max lifespan", {
  
  h <- Household$new()

  expect_error(
    h$get_lifespan(current_date = "2020-07-15")
  )
  
  h$add_member(
    HouseholdMember$new(
      name       = "test_member_older",  
      birth_date = "1980-07-15"
    )  
  )  
  expect_equal(
    h$get_lifespan(current_date = "2020-07-15"),
    60
  )

  h$add_member(
    HouseholdMember$new(
      name       = "test_member_younger",  
      birth_date = "1990-07-15"
    )
  )  
  expect_equal(
    h$get_lifespan(current_date = "2020-07-15"),
    70
  )  

  h$set_lifespan(90)
  expect_equal(
    h$get_lifespan(current_date = "2020-07-15"),
    90
  )
})

test_that("setting planned income", {

  household <- Household$new()
  household$add_member(
    HouseholdMember$new(
      name       = "test_member_older",  
      birth_date = "1980-07-15"
    )  
  )  
  household$add_member(
    HouseholdMember$new(
      name       = "younger",  
      birth_date = "1990-07-15"
    )
  )  
  
  expect_true(is.list(household$expected_income))
  expect_equal(
    length(household$expected_income),
    0
  )

  household$expected_income <- list(
    "income1" = c(
      "members$older$age >= 44 & members$older$age < 46 ~ 100",
      "members$older$age >= 46 ~ 300"
    ),
    "income2" = c(
      "members$younger$age >= 34 & members$younger$age < 36 ~ 44",
      "members$younger$age >= 36 ~ 55"
    )
  )
  expect_equal(
    length(household$expected_income),
    2
  )
})

test_that("setting planned non-discretionary spending", {

  household <- Household$new()
  household$add_member(
    HouseholdMember$new(
      name       = "test_member_older",  
      birth_date = "1980-07-15"
    )  
  )  
  household$add_member(
    HouseholdMember$new(
      name       = "younger",  
      birth_date = "1990-07-15"
    )
  )  
  
  expect_true(is.list(household$expected_spending))
  expect_equal(
    length(household$expected_spending),
    0
  )

  household$expected_spending <- list(
    "spending1" = c(
      "members$older$age >= 44 & members$older$age < 46 ~ 100",
      "members$older$age >= 46 ~ 300"
    ),
    "spending2" = c(
      "members$younger$age >= 34 & members$younger$age < 36 ~ 44",
      "members$younger$age >= 36 ~ 55"
    )
  )
  expect_equal(
    length(household$expected_spending),
    2
  )
})

test_that("setting household risk tolerance", {

  household <- Household$new()
  expect_equal(household$risk_tolerance, 0.5)

  household$risk_tolerance <- 0.35
  expect_equal(household$risk_tolerance, 0.35)
})

test_that("setting household consumption impatience preference", {

  household <- Household$new()
  expect_equal(household$consumption_impatience_preference, 0.04)

  household$consumption_impatience_preference <- 0.08
  expect_equal(household$consumption_impatience_preference, 0.08)
})

test_that("setting household smooth consumption preference", {

  household <- Household$new()
  expect_equal(household$smooth_consumption_preference, 1)

  household$smooth_consumption_preference <- 0.5
  expect_equal(household$smooth_consumption_preference, 0.5)
})

test_that("calculating joint Gompertz parameters for 1 member", {

  test_birth_date   <- "1955-07-15"
  test_current_date <- "2020-07-15"
  
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  members$mode       <- 80
  members$dispersion <- 10

  household <- Household$new()
  household$add_member(members)
  
  survival <- household$calc_survival(current_date = test_current_date)
  expect_equal(survival$mode, 80, tolerance = 0.5)
  expect_equal(survival$dispersion, 10, tolerance = 0.5)
  expect_equal(
    survival$data |> 
    dplyr::filter(year == 85 - 65) |> 
    dplyr::pull(gompertz),
    0.2396, 
    tolerance = 0.001
  )
})

test_that("calculating joint Gompertz parameters for 2 members", {
  
  test_current_date <- "2020-01-01"
  hm1 <- 
    HouseholdMember$new(
      name       = "member1",
      birth_date = "1955-01-01"
    )
  expect_equal(
    hm1$calc_age(current_date = test_current_date), 
    65, 
    tolerance = 0.01
  )
  hm1$mode       <- 88
  hm1$dispersion <- 10.65

  hm2 <- 
    HouseholdMember$new(
      name = "member2",
      birth_date = "1955-01-01"
    )
  expect_equal(
    hm2$calc_age(current_date = test_current_date), 
    65, 
    tolerance = 0.01
  )
  hm2$mode       <- 91
  hm2$dispersion <- 8.88

  household <- Household$new()
  household$add_member(hm1)
  household$add_member(hm2)
  household$set_lifespan(45)
  
  params <- household$calc_survival(current_date = test_current_date) 
expect_equal(
     params$mode, 
    93.22423,
    tolerance = 1e-5
   )
   expect_equal(
     params$dispersion, 
    5.38006535,
    tolerance = 1e-8
   )
})

test_that("calculating joint Gompertz parameters for 3 members", {
  
  test_current_date <- "2020-01-01"
  hm1 <- 
    HouseholdMember$new(
      name       = "member1",
      birth_date = "1955-01-01"
    )
  expect_equal(
    hm1$calc_age(current_date = test_current_date), 
    65, 
    tolerance = 0.01
  )
  hm1$mode       <- 88
  hm1$dispersion <- 10.65

  hm2 <- 
    HouseholdMember$new(
      name = "member2",
      birth_date = "1955-01-01"
    )
  expect_equal(
    hm2$calc_age(current_date = test_current_date), 
    65, 
    tolerance = 0.01
  )
  hm2$mode       <- 91
  hm2$dispersion <- 8.88

  hm3 <- 
    HouseholdMember$new(
      name = "member3",
      birth_date = "1955-01-01"
    )
  expect_equal(
    hm3$calc_age(current_date = test_current_date), 
    65, 
    tolerance = 0.01
  )
  hm3$mode       <- 95
  hm3$dispersion <- 8.88

  household <- Household$new()
  household$add_member(hm1)
  household$add_member(hm2)
  household$add_member(hm3)
  household$set_lifespan(45)
 
  params <- household$calc_survival(current_date = test_current_date)

  expect_true(is.double(params$data[[hm1$get_name()]]))
  expect_true(is.double(params$data[[hm2$get_name()]]))
  expect_true(is.double(params$data[[hm3$get_name()]]))

  expect_equal(
    params$mode, 
    95.691272,
    tolerance = 1e-6
  )
  expect_equal(
    params$dispersion, 
    3.8796684,
    tolerance = 1e-7
  )
})

test_that("getting min_age - age of the youngest member", {

  test_current_date <- "2020-01-01"
  hm1 <- 
    HouseholdMember$new(
      name       = "member1",
      birth_date = "1955-01-01"
    )
  hm2 <- 
    HouseholdMember$new(
      name = "member2",
      birth_date = "1965-01-01"
    )
  hm3 <- 
    HouseholdMember$new(
      name = "member3",
      birth_date = "1975-01-01"
    )
  household <- Household$new()
  household$add_member(hm1)
  household$add_member(hm2)
  household$add_member(hm3) 

  hm1$calc_age(current_date = test_current_date)
  hm2$calc_age(current_date = test_current_date)
  hm3$calc_age(current_date = test_current_date)

  expect_equal(
    household$get_min_age(current_date = test_current_date),
    45,
    tolerance = 0.01
  )
})

test_that("cloning works", {

  skip_if_not(interactive())
  skip_on_ci()
  skip_on_cran()

  test_birth_date   <- "1955-07-15"
  members <- HouseholdMember$new(
    name       = "test_name",
    birth_date = test_birth_date
  )
  members$set_event("retirement", 65)
  household <- Household$new()
  household$add_member(members)
  expect_equal(
    household$get_members()$test_name$get_events()$retirement$start_age,
    65
  )

  household_hash <- rlang::hash(household)
  expect_equal(
    rlang::hash(household), 
    household_hash
  )

  household_bis <- unserialize(serialize(household, NULL))

  expect_equal(
    rlang::hash(household_bis), 
    household_hash
  )
  expect_equal(
    rlang::hash(household), 
    household_hash
  )

  household_bis$get_members()$test_name$set_event("retirement", 100)
  expect_equal(
    household_bis$get_members()$test_name$get_events()$retirement$start_age,
    100
  )
  expect_equal(
    household$get_members()$test_name$get_events()$retirement$start_age,
    65
  )

  expect_equal(rlang::hash(household), household_hash)
})
