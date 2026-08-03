calc_present_value <- function(
  cashflow,
  discount_rate
) {
  v <- time_value_discount <- NULL

  n_years        <- length(cashflow)
  present_values <- vector("numeric", n_years)

  # Create a matrix of discount factors
  discount_matrix <- matrix(0, nrow = n_years, ncol = n_years)
  for (i in 1:n_years) {
    for (j in i:n_years) {
      discount_matrix[i, j] <- 1 / ((1 + discount_rate)^(j - i))
    }
  }

  # Create a sparse matrix representation with just the relevant values
  cashflow_matrix <- matrix(0, nrow = n_years, ncol = n_years)
  for (i in 1:n_years) {
    cashflow_matrix[i, i:n_years] <- cashflow[i:n_years]
  }

  # Calculate present values as the sum of discounted cashflows
  present_values <- rowSums(discount_matrix * cashflow_matrix)
  present_values
}
