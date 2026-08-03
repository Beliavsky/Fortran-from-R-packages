balsh <- function(FA, CA, INV, FL, CL) {
  TA <- FA + CA # Total Assets
  TL <- FL + CL # Total Liabilities
  Leverage <- TL / TA # Leverage (L/A)
  EQ <- TA - TL # Equity
  DEratio <- TL / EQ # Debt to Equity Ratio
  WCAP <- CA - CL # Working Capital
  CURatio <- CA / CL # Current Ratio
  QUICK <- (CA - INV) / CL # Quick (or acid-test Ratio)
  
  results <- data.frame(
    Metric = c("Working Capital", "Current Ratio", "Acid Test Ratio", "Leverage (L/A) Ratio", "Debt to Equity (D/E) Ratio"),
    Value = sprintf("%.2f", c(WCAP, CURatio, QUICK, Leverage, DEratio))
  )
  
  return(results)
}



