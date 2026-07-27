## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(greeks)

## -----------------------------------------------------------------------------
european_put <- BS_European_Greeks(
  initial_price = 120,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 4.5,
  dividend_yield = 0.015,
  volatility = 0.22,
  payoff = "put",
  greek = c("fair_value", "delta", "gamma", "vega", "theta", "rho")
)

round(european_put, 4)

## -----------------------------------------------------------------------------
round(
  Greeks(
    initial_price = 120,
    exercise_price = 100,
    r = 0.02,
    time_to_maturity = 4.5,
    dividend_yield = 0.015,
    volatility = 0.22,
    payoff = "put",
    greek = c("fair_value", "delta", "gamma")
  ),
  4
)

## -----------------------------------------------------------------------------
digital_call <- BS_European_Greeks(
  initial_price = 100,
  exercise_price = 105,
  r = 0.03,
  time_to_maturity = 1,
  volatility = 0.25,
  payoff = "cash_or_nothing_call",
  greek = c("fair_value", "delta", "vega")
)

round(digital_call, 4)

## -----------------------------------------------------------------------------
base_args <- list(
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1.5,
  dividend_yield = 0,
  volatility = 0.3,
  payoff = "call"
)

fair_value_at <- function(initial_price) {
  do.call(
    BS_European_Greeks,
    c(base_args, list(initial_price = initial_price, greek = "fair_value"))
  )
}

step_size <- 1e-4
finite_difference_delta <-
  (fair_value_at(100 + step_size) - fair_value_at(100 - step_size)) /
  (2 * step_size)

exact_delta <- do.call(
  BS_European_Greeks,
  c(base_args, list(initial_price = 100, greek = "delta"))
)

round(
  c(
    exact_delta = exact_delta,
    finite_difference_delta = finite_difference_delta,
    absolute_error = abs(exact_delta - finite_difference_delta)
  ),
  8
)

## -----------------------------------------------------------------------------
greeks_to_compare <- c("fair_value", "delta", "vega", "theta", "rho", "gamma")

exact <- BS_European_Greeks(
  initial_price = 110,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  volatility = 0.25,
  payoff = "call",
  greek = greeks_to_compare
)

monte_carlo <- Malliavin_European_Greeks(
  initial_price = 110,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  volatility = 0.25,
  payoff = "call",
  greek = greeks_to_compare,
  paths = 50000,
  seed = 42,
  antithetic = TRUE
)

round(rbind(exact = exact, malliavin_monte_carlo = monte_carlo), 4)

