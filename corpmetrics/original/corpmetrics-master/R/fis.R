fis <- function(FV, CR, YTM, MAT, SEMI = FALSE) {
  # Checking if payments are semi-annual (default is annual payments)
  if (SEMI) {
    n <- 2 * MAT # n is periods
    c <- CR / 2 # c is periodic coupon rate
    y <- YTM / 2 # y is semi-annual Yield to Maturity
  } else {
    n <- MAT # Maturity
    c <- CR # Coupon Rate
    y <- YTM # Yield to Maturity
  }

  # Calculate the price
  price <- 0
  for (t in 1:n) {
    price <- price + (FV * c) / ((1 + y) ^ t) # Discounted Cash Flows
  }
  price <- price + FV / ((1 + y) ^ n) # Add the Discounted face value payment at maturity

  # Calculate the Macaulay Duration
  duration <- 0
  for (t in 1:n) {
    duration <- duration + (t * (FV * c)) / ((1 + y) ^ t)
  }
  duration <- duration + (n * FV) / ((1 + y) ^ n)
  duration <- duration / price

  if (SEMI) {
    duration <- duration / 2
  }

  modified_duration <- duration / (1 + YTM)

  results_df <- data.frame(
    Description = c("Price", "Macaulay Duration", "Modified Duration"),
    Value = c(
      sprintf("%.2f Currency Units", price),
      sprintf("%.2f Years", duration),
      sprintf("%.2f Years", modified_duration)
    )
  )

  return(results_df)
}



