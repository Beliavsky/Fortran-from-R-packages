# Run help(DNF) for description and argument details.
DNF <- function(dynamics, ...) UseMethod("DNF") 

DNF.default <- function(dynamics, ...){
  stop("This class of state-space models is not supported by the DNF function.")
}

DNF.dynamicsSVM <- function(dynamics, data, factors = NULL, N = 50, K = 20, R = 1, grids = "Default", ...) {
  # If explanatory variables or regression coefficients were passed,
  # check matching dates, dimensions, etc...
  if(!is.null(factors) | !is.null(dynamics$coefs)){
    factorsModelsCheck(coefs = dynamics$coefs,
                       factors = factors, data = data)
  }
  
  if(is.xts(data)){
    # Extract returns to run DNF with
    ret <- coredata(data)
    # Transform column vector of returns into row vectors
    if(nrow(ret) > 1){
      ret <- t(ret)
      if(nrow(ret) > 1){
        stop("You should pass a vector of returns (i.e., 1xT or Tx1 dimensional vector).")
      }
    }
  }else{
    ret <- data
    
  }
  T <- length(ret)
  
  if(!is.null(factors)){
    if(T != nrow(factors)){
      # Check that factors has the correct dimensions.
      stop("Your factor matrix has incorrect dimensions. Verify that either its number of rows match the length or your returns data.")
    }
    
    ret <- ret - t(factors %*% as.matrix(dynamics$coefs))
    }

  likelihoods <- seq(from = 0, to = 0, length = T)
  log_likelihood <- 0
  # Creating grid for listed models
  if (length(grids) == 1) { # Check if users input a custom grid
    grids <- gridMaker(dynamics, R = R, N = N, K = K)
  }
  # Get the lengths from the grids list (as the it could be provided by the user)
  var_mid_points <- grids$var_mid_points
  N <- length(var_mid_points)
  R <- length(grids$j_nums)
  K <- length(grids$jump_mid_points)

  # Set a matrix to store the filtering distribution at each time step
  filter_grid <- matrix(0, nrow = N, ncol = T + 1)
  filter_grid[(N:1), 1] <- rep(1 / N, times = N) # Uniform prior

  d_probs <- probCalculator(grids = grids, R = R, N = N, K = K, data = ret, dynamics = dynamics)
    
    # Function to compute measurement equation times the other probabilities.
  for (t in (1:T)) {
    # Time t Likelihood contribution:
    probs <- Cpp_prodfun(d_probs[, t], filter_grid[(N:1), t]) # C++ function to multiply the previous filtering distribution with the appropriate probabilities.
    likelihoods[t] <- sum(probs)
    l = likelihoods[t]
    
    # Get the posterior filtering distribution:
    if(is.na(l) | is.nan(l)){
      return(list(log_likelihood = -Inf, filter_grid = filter_grid,
                  likelihoods = likelihoods, grids = grids, dynamics = dynamics))
    }
    if(l > 0){
      filter_grid[(N:1), t + 1] <- Cpp_rowSums_modN(probs, N) / l
          }
  else{
    filter_grid[, t + 1] <- filter_grid[, t]
  }
  }
  # Compute the likelihood:
  log_likelihood <- sum(log(likelihoods))

  SVDNF <- list(log_likelihood = log_likelihood, filter_grid = filter_grid,
                likelihoods = likelihoods, grids = grids, dynamics = dynamics, data = data)
  class(SVDNF) = "SVDNF"
  return(SVDNF)
}

# function to ensure that factor models passes a series of checks
factorsModelsCheck <- function(coefs, factors, data){
  # Check if both coefficients and factors have been passed.
  if(!is.null(coefs) & is.null(factors)){
    stop("You must pass series  of factors/explanatory variables as there are regression coefficients in the model dynamics. If you do not want explanatory variables, remove the coefs from the model dynamics.")
  }
  if(is.null(coefs) & !is.null(factors)){
    stop("You must have coefficients in your model dynamics as there explanatory variables passed to the DNF function. If you do not want explanatory variables, keep factors = NULL in the DNF function.")
  }
  
  # if factors is also an xts...
  if(is.xts(factors) & is.xts(data)){
    # check that the dates match
    if(!(index(factors) == index(data))){
      stop("The dates for the explanatory variables (factors) and for the returns (data) do not match.")
    }
  }
  # Number of coefs matches number of explanatory
  if(length(coefs) != ncol(factors)){
    stop("The number of regression coefficients does not match the number of factors passed the DNF function (check ncol(factors) and the length of your coefficient vectors, coefs in the dynamics object.)")
  }
  
  # Length of returns series matches length of explanatory variables.
  if(length(data)!= nrow(factors)){
    stop("The length of the returns series does not match the number of rows in the factors matrix.")
  }
}

print.SVDNF <- function(x, ...){
  cat("\nModel:\n")
  cat(x$dynamics$model, '\n')
  # MLE coefficient estimates
  cat("\nLog-likelihood:\n")
  print(x$log_likelihood)
  cat("\n")
  invisible(x)

}
