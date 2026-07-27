## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(greeks)

## -----------------------------------------------------------------------------
true_volatility <- 0.28

option_price <- BS_European_Greeks(
  initial_price = 100,
  exercise_price = 105,
  r = 0.03,
  time_to_maturity = 1.25,
  dividend_yield = 0.01,
  volatility = true_volatility,
  payoff = "call",
  greek = "fair_value"
)

implied_volatility <- BS_Implied_Volatility(
  option_price = option_price,
  initial_price = 100,
  exercise_price = 105,
  r = 0.03,
  time_to_maturity = 1.25,
  dividend_yield = 0.01,
  payoff = "call",
  start_volatility = 0.2
)

round(
  c(
    option_price = option_price,
    true_volatility = true_volatility,
    implied_volatility = implied_volatility
  ),
  6
)

## -----------------------------------------------------------------------------
reconstructed_price <- BS_European_Greeks(
  initial_price = 100,
  exercise_price = 105,
  r = 0.03,
  time_to_maturity = 1.25,
  dividend_yield = 0.01,
  volatility = implied_volatility,
  payoff = "call",
  greek = "fair_value"
)

round(
  c(
    original_price = option_price,
    reconstructed_price = reconstructed_price,
    absolute_error = abs(option_price - reconstructed_price)
  ),
  10
)

## -----------------------------------------------------------------------------
geometric_true_volatility <- 0.35

geometric_price <- Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = geometric_true_volatility,
  option_type = "Geometric Asian",
  payoff = "put",
  greek = "fair_value"
)

geometric_implied_volatility <- Implied_Volatility(
  option_price = geometric_price,
  initial_price = 100,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  dividend_yield = 0,
  option_type = "Geometric Asian",
  payoff = "put"
)

geometric_reconstructed_price <- Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = geometric_implied_volatility,
  option_type = "Geometric Asian",
  payoff = "put",
  greek = "fair_value"
)

round(
  c(
    original_price = geometric_price,
    true_volatility = geometric_true_volatility,
    implied_volatility = geometric_implied_volatility,
    reconstructed_price = geometric_reconstructed_price,
    absolute_error = abs(geometric_price - geometric_reconstructed_price)
  ),
  6
)

## -----------------------------------------------------------------------------
near_zero_volatility_price <- Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.02,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = 1e-12,
  option_type = "Geometric Asian",
  payoff = "put",
  greek = "fair_value"
)

round(near_zero_volatility_price, 6)

