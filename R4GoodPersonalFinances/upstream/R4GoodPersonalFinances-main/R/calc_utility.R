calc_utility <- function(x, parameter) {

  x <- ifelse(x < 0, 0, x)

  utility <- numeric(length(x))
  
  utility[parameter == 1] <- log(x[parameter == 1])

  utility[parameter != 1] <- 
    (parameter[parameter != 1] / (parameter[parameter != 1] - 1)) * 
    (
      x[parameter != 1] ^ (
        (parameter[parameter != 1] - 1) / parameter[parameter != 1]
      ) - 1
    )
  
  utility[utility == -Inf] <- 0
  utility[is.nan(utility)] <- 0
  utility
}

calc_inverse_utility <- function(
  utility,
  parameter
) {
  
  inverse_utility <- numeric(length(utility))
  
  inverse_utility[parameter == 1] <- exp(utility[parameter == 1])

  inverse_utility[parameter != 1] <- 
    (
      ((parameter[parameter != 1] - 1) / parameter[parameter != 1]) * 
      utility[parameter != 1] + 1
    ) ^ (
      parameter[parameter != 1] / (parameter[parameter != 1] - 1)
    )
  
  inverse_utility
}
