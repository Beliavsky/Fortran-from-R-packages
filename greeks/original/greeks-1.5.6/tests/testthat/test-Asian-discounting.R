test_that("geometric Asian prices use the risk-free discount rate", {
  initial_price <- 100
  exercise_price <- 100
  r <- 0.05
  dividend_yield <- 0.03
  time_to_maturity <- 2
  volatility <- 0.2

  log_mean <-
    log(initial_price) +
    0.5 * (r - dividend_yield - volatility^2 / 2) * time_to_maturity
  log_variance <- volatility^2 * time_to_maturity / 3
  log_sd <- sqrt(log_variance)
  d_put <- (log_mean - log(exercise_price)) / log_sd
  d_call <- d_put + log_sd
  discount_factor <- exp(-r * time_to_maturity)

  expected_call <-
    discount_factor *
    (exp(log_mean + log_variance / 2) * pnorm(d_call) -
       exercise_price * pnorm(d_put))
  expected_put <-
    discount_factor *
    (exercise_price * pnorm(-d_put) -
       exp(log_mean + log_variance / 2) * pnorm(-d_call))

  actual_call <- BS_Geometric_Asian_Greeks(
    initial_price = initial_price,
    exercise_price = exercise_price,
    r = r,
    time_to_maturity = time_to_maturity,
    volatility = volatility,
    dividend_yield = dividend_yield,
    payoff = "call",
    greek = "fair_value"
  )
  actual_put <- BS_Geometric_Asian_Greeks(
    initial_price = initial_price,
    exercise_price = exercise_price,
    r = r,
    time_to_maturity = time_to_maturity,
    volatility = volatility,
    dividend_yield = dividend_yield,
    payoff = "put",
    greek = "fair_value"
  )

  expect_equal(unname(actual_call), expected_call, tolerance = 1e-12)
  expect_equal(unname(actual_put), expected_put, tolerance = 1e-12)
})

test_that("Asian estimators discount independently of the cost of carry", {
  time_to_maturity <- 0.75
  rate_shift <- 0.03125
  discount_scale <- exp(-rate_shift * time_to_maturity)

  common_args <- list(
    initial_price = 100,
    exercise_price = 100,
    r = 0.0625,
    time_to_maturity = time_to_maturity,
    volatility = 0.25,
    dividend_yield = 0.03125,
    payoff = "call",
    steps = 8,
    paths = 256,
    seed = 42
  )

  shifted_args <- common_args
  shifted_args$r <- shifted_args$r + rate_shift
  shifted_args$dividend_yield <-
    shifted_args$dividend_yield + rate_shift

  expect_discount_shift <- function(fun, greek) {
    base_args <- common_args
    base_args$greek <- greek
    comparison_args <- shifted_args
    comparison_args$greek <- greek

    base_value <- do.call(fun, base_args)
    shifted_value <- do.call(fun, comparison_args)
    expected <- discount_scale * base_value

    theta_names <- intersect(c("theta", "theta_d"), names(base_value))
    expected[theta_names] <-
      discount_scale *
      (base_value[theta_names] + rate_shift * base_value["fair_value"])

    expect_equal(shifted_value, expected, tolerance = 1e-10)
  }

  expect_discount_shift(
    BS_Malliavin_Asian_Greeks,
    c("fair_value", "delta", "rho", "vega")
  )
  expect_discount_shift(
    Malliavin_Asian_Greeks,
    c("fair_value", "theta", "theta_d")
  )
  expect_discount_shift(
    Malliavin_Geometric_Asian_Greeks,
    c("fair_value", "delta", "rho", "vega", "theta", "gamma")
  )
})
