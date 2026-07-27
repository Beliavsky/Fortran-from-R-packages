## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(greeks)

## -----------------------------------------------------------------------------
american_put <- Binomial_American_Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.03,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = 0.3,
  payoff = "put",
  greek = c("fair_value", "delta", "gamma", "vega", "theta", "rho"),
  steps = 200
)

round(american_put, 4)

## -----------------------------------------------------------------------------
round(
  Greeks(
    initial_price = 100,
    exercise_price = 100,
    r = 0.03,
    time_to_maturity = 1,
    dividend_yield = 0,
    volatility = 0.3,
    payoff = "put",
    option_type = "American",
    greek = c("fair_value", "delta", "rho"),
    steps = 200
  ),
  4
)

## -----------------------------------------------------------------------------
european_put_value <- BS_European_Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.03,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = 0.3,
  payoff = "put",
  greek = "fair_value"
)

american_put_value <- Binomial_American_Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.03,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = 0.3,
  payoff = "put",
  greek = "fair_value",
  steps = 400
)

round(
  c(
    european_put = european_put_value,
    american_put = american_put_value,
    early_exercise_premium = american_put_value - european_put_value
  ),
  4
)

## -----------------------------------------------------------------------------
steps <- c(25, 50, 100, 200, 400)

prices <- vapply(
  steps,
  function(n_steps) {
    Binomial_American_Greeks(
      initial_price = 100,
      exercise_price = 100,
      r = 0.03,
      time_to_maturity = 1,
      dividend_yield = 0,
      volatility = 0.3,
      payoff = "put",
      greek = "fair_value",
      steps = n_steps
    )
  },
  numeric(1)
)

data.frame(steps = steps, fair_value = round(prices, 4))

## -----------------------------------------------------------------------------
fair_value_at <- function(initial_price) {
  Binomial_American_Greeks(
    initial_price = initial_price,
    exercise_price = 100,
    r = 0.03,
    time_to_maturity = 1,
    dividend_yield = 0,
    volatility = 0.3,
    payoff = "put",
    greek = "fair_value",
    steps = 200
  )
}

step_size <- 0.01
finite_difference_delta <-
  (fair_value_at(100 + step_size) - fair_value_at(100 - step_size)) /
  (2 * step_size)

reported_delta <- Binomial_American_Greeks(
  initial_price = 100,
  exercise_price = 100,
  r = 0.03,
  time_to_maturity = 1,
  dividend_yield = 0,
  volatility = 0.3,
  payoff = "put",
  greek = "delta",
  steps = 200,
  eps = step_size
)

round(
  c(
    reported_delta = reported_delta,
    finite_difference_delta = finite_difference_delta,
    absolute_error = abs(reported_delta - finite_difference_delta)
  ),
  8
)

