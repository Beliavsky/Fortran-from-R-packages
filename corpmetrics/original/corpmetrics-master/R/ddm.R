ddm <- function(DIV, RETURN, G1 = NULL, G2 = NULL, PER = NULL) {
  if (is.null(G1) && is.null(G2) && is.null(PER)) {
    # Zero Growth Model
    value <- DIV / RETURN
    method <- "Zero Growth Model"
    
  } else if (!is.null(G1) && is.null(G2) && is.null(PER)) {
    # Gordon's Model
    div1 <- DIV * (1 + G1)
    value <- div1 / (RETURN - G1)
    method <- "Gordon's Model"
    
  } else if (!is.null(G1) && !is.null(G2) && !is.null(PER)) {
    # Differential Growth Model
    div1 <- DIV * (1 + G1)
    value1 <- div1 / (RETURN - G1) *
      (1 - (1 + G1) ^ PER / (1 + RETURN) ^ PER)
    div_n1 <- DIV * (1 + G1) ^ (PER + 1)
    value2 <- div_n1 / (RETURN - G2) / (1 + RETURN) ^ PER
    value <- value1 + value2
    method <- "Differential Growth Model"
    
  } else {
    stop("Invalid combination of parameters. Please provide the appropriate inputs.")
  }
  
  results_df <- data.frame(
    Description = c("Model Employed", "Stock's Value"),
    Value = c(method, sprintf("%.2f Currency Units", value))
  )
  
  return(results_df)
}



