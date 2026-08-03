loan <- function(AMOUNT, RATE, PER) {
  installment <- AMOUNT * RATE / (1 - (1 + RATE)^(-PER))
  loan_balance <- AMOUNT
  table <- data.frame(
    Period = 1:PER,
    Installment = rep(installment, PER),
    Interest = numeric(PER),
    Principal = numeric(PER),
    Balance = numeric(PER)
  )
  
  for (period in 1:PER) {
    payment <- loan_balance * RATE
    principal_payment <- installment - payment
    loan_balance <- loan_balance - principal_payment
    
    table[period, c("Interest", "Principal", "Balance")] <- c(
      round(payment, 2),
      round(principal_payment, 2),
      round(loan_balance, 2)
    )
  }
  
  # Format numeric columns to remove scientific notation
  table$Installment <- format(round(table$Installment, 2), nsmall = 2)
  table$Interest <- format(round(table$Interest, 2), nsmall = 2)
  table$Principal <- format(round(table$Principal, 2), nsmall = 2)
  table$Balance <- format(round(table$Balance, 2), nsmall = 2)
  
  # Results data frame
  results_df <- data.frame(
    Description = c("Loan Installment", "Total Repayment"),
    Value = c(
      sprintf("%.2f Currency Units", installment),
      sprintf("%.2f Currency Units", installment * PER)
    )
  )
  
  return(list(Summary = results_df, AmortizationTable = table))
}




