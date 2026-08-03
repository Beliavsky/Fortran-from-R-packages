insta <- function(REV, COS, NET, PREF = NULL, SHARES = NULL, PPS = NULL) {
  gross_profit_margin <- ((REV - COS) / REV) * 100
  net_profit_margin <- (NET / REV) * 100
  if (!is.null(PREF) && !is.null(SHARES)) {
    eps <- (NET - PREF) / SHARES
  } else {
    eps <- NA
  }
  if (!is.null(eps) && !is.null(PPS)) {
    pe_ratio <- PPS / eps
  } else {
    pe_ratio <- NA
  }
  
  results_df <- data.frame(
    Metric = c("Gross Profit Margin (%)", "Net Profit Margin (%)", "Earnings Per Share (EPS)", "Price to Earnings (P/E) Ratio"),
    Value = c(
      sprintf("%.2f", gross_profit_margin),
      sprintf("%.2f", net_profit_margin),
      ifelse(!is.na(eps), sprintf("%.2f", eps), NA),
      ifelse(!is.na(pe_ratio), sprintf("%.2f", pe_ratio), NA)
    )
  )
  
  return(results_df)
}


