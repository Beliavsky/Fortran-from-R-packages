# Shared test fixtures

# Realistic US Treasury yields (normal curve)
make_us_yields <- function() {
  list(
    maturities = c(0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 30),
    rates = c(0.052, 0.050, 0.048, 0.045, 0.043, 0.042, 0.041,
              0.040, 0.042, 0.043)
  )
}

# Flat yield curve
make_flat_yields <- function() {
  list(
    maturities = c(1, 2, 5, 10, 30),
    rates = rep(0.04, 5)
  )
}

# Inverted yield curve
make_inverted_yields <- function() {
  list(
    maturities = c(0.25, 0.5, 1, 2, 5, 10, 30),
    rates = c(0.055, 0.054, 0.052, 0.048, 0.043, 0.040, 0.038)
  )
}

# Simple 4-point curve for quick tests
make_simple_yields <- function() {
  list(
    maturities = c(2, 5, 10, 30),
    rates = c(0.045, 0.042, 0.040, 0.043)
  )
}
