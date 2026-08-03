capm <- function(Rf, Ri, Rm) {
  # Calculate expected returns
  expected_asset <- mean(Ri) # Expected asset's returns
  expected_market <- mean(Rm) # Expected market's returns
  riskfree <- mean(Rf) # Risk-free rate

  # Calculate beta
  beta <- cov(Ri, Rm) / var(Rm)

  # Calculate required return based on CAPM
  required_return <- riskfree + beta * (expected_market - riskfree)

  # Determine valuation status
  if (required_return < expected_asset) {
    value <- "undervalued"
  } else if (required_return > expected_asset) {
    value <- "overvalued"
  } else {
    value <- "fairly valued"
  }

  # results data frame
  results_df <- data.frame(
    Description = c("Required Return based on CAPM", "Expected Return of the Security", "Stock's Beta", "Valuation Status"),
    Value = c(
      sprintf("%.2f", required_return),
      sprintf("%.2f", expected_asset),
      sprintf("%.2f", beta),
      value
    )
  )

  return(results_df)
}


