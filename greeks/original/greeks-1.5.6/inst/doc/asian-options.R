## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(greeks)

## -----------------------------------------------------------------------------
geometric_exact <- BS_Geometric_Asian_Greeks(
  initial_price = 110,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1.5,
  dividend_yield = 0.01,
  volatility = 0.25,
  payoff = "call",
  greek = c("fair_value", "delta", "rho", "vega", "theta", "gamma")
)

round(geometric_exact, 4)

## -----------------------------------------------------------------------------
geometric_price_at <- function(volatility) {
  BS_Geometric_Asian_Greeks(
    initial_price = 110,
    exercise_price = 100,
    r = 0.02,
    time_to_maturity = 1.5,
    dividend_yield = 0.01,
    volatility = volatility,
    payoff = "call",
    greek = "fair_value"
  )
}

step_size <- 1e-4
finite_difference_vega <-
  (geometric_price_at(0.25 + step_size) -
     geometric_price_at(0.25 - step_size)) /
  (2 * step_size)

round(
  c(
    exact_vega = geometric_exact["vega"],
    finite_difference_vega = finite_difference_vega,
    absolute_error = abs(geometric_exact["vega"] - finite_difference_vega)
  ),
  8
)

## -----------------------------------------------------------------------------
greeks_to_compare <- c("fair_value", "delta", "rho", "vega")

geometric_monte_carlo <- Malliavin_Geometric_Asian_Greeks(
  initial_price = 110,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1.5,
  dividend_yield = 0.01,
  volatility = 0.25,
  payoff = "call",
  greek = greeks_to_compare,
  paths = 20000,
  steps = 24,
  seed = 42,
  antithetic = TRUE
)

round(
  rbind(
    exact = geometric_exact[greeks_to_compare],
    malliavin_monte_carlo = geometric_monte_carlo
  ),
  4
)

## -----------------------------------------------------------------------------
arithmetic_asian <- BS_Malliavin_Asian_Greeks(
  initial_price = 110,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1.5,
  dividend_yield = 0.01,
  volatility = 0.25,
  payoff = "call",
  greek = c("fair_value", "delta", "rho", "vega"),
  paths = 10000,
  steps = 24,
  seed = 42
)

round(arithmetic_asian, 4)

## -----------------------------------------------------------------------------
round(
  Greeks(
    initial_price = 110,
    exercise_price = 100,
    r = 0.02,
    time_to_maturity = 1.5,
    dividend_yield = 0.01,
    volatility = 0.25,
    payoff = "call",
    option_type = "Asian",
    greek = c("fair_value", "delta", "rho", "vega")
  ),
  4
)

## -----------------------------------------------------------------------------
price_grid <- Malliavin_Geometric_Asian_Greeks(
  initial_price = c(95, 100, 105),
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  volatility = 0.25,
  payoff = "put",
  greek = c("fair_value", "delta"),
  paths = 5000,
  steps = 12,
  seed = 42
)

round(price_grid, 4)

## -----------------------------------------------------------------------------
digital_call <- function(x, exercise_price) {
  ifelse(x >= exercise_price, 1, 0)
}

invisible(capture.output(
  custom_digital <- Malliavin_Geometric_Asian_Greeks(
    initial_price = 100,
    exercise_price = 100,
    r = 0.02,
    time_to_maturity = 1,
    volatility = 0.25,
    payoff = digital_call,
    greek = c("fair_value", "delta"),
    paths = 5000,
    steps = 12,
    seed = 42
  )
))

round(custom_digital, 4)

