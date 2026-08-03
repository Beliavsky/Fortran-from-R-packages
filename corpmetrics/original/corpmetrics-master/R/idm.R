idm <- function(CFS, COST) {
  if (length(CFS) != length(COST)) {
    stop("The length of Cash Flows and Cost of Capital must be the same")
  }
  
  npv <- function(CFS, COST) {
    sum(CFS / (1 + COST)^(0:(length(CFS) - 1)))
  }
  
  npv_result <- npv(CFS, COST)
  npv_function <- function(COST) npv(CFS, COST)
  irr_result <- uniroot(npv_function, c(-1, 1))$root
  irr_result <- irr_result * 100
  
  results_df <- data.frame(
    Method = c("Net Present Value (NPV)", "Internal Rate of Return (IRR)"),
    Value = c(
      sprintf("%.2f Currency Units", npv_result),
      sprintf("%.2f %%", irr_result)
    )
  )
  
  return(results_df)
}






