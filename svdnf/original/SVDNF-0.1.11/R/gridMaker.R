# Function gridMaker:
#   This function creates grids for the pre-specified models in the SVDNF package.
#
# Inputs:
#    - N, K, R : Sizes for the volatility, volatility jumps, and number of jumps grids.
#   	- models parameters : List of parameters used in the
#                          pre-specified models (see help(DNF) for a description).
#    - model : Model for which we are applying the DNF.
# Output:
#   	- grids : A list that contains the sequence of volatility nodes, var_mid_points,
#              the nodes for the number of jumps, j_nums, and
#              the volatility jump size nodes, jump_mid_points.

gridMaker <- function(model_dynamics, ...) UseMethod("gridMaker") 

gridMaker.default <- function(model_dynamics, ...){
  stop("Please select a built-in model or input your own grid to the DNF")
}
# Model from Duffie, Pan, and Singleton (2000):
gridMaker.DuffiePanSingleton <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  nu <- model_dynamics$nu
  sigma <- unlist(sigma_x_params[1])
  theta <- unlist(mu_x_params[2]); kappa <- unlist(mu_x_params[1])
  var_mid_points <- seq(from = sqrt(max(theta - (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa), 0.0000001)), to = max(sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), sqrt(0.05)), length = N)^2
  j_nums <- seq(from = 0, to = R, by = 1)
  jump_mid_points <- seq(from = 0.000001, to = (3 + log(K)) * sqrt(R) * nu, length = K)
  
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}
# Model from Bates (1996):
gridMaker.Bates <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  sigma <- unlist(sigma_x_params[1])
  theta <- unlist(mu_x_params[2]); kappa <- unlist(mu_x_params[1])
  
  var_mid_points <- seq(from = sqrt(max(theta - (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa), 0.0000001)), to = max(sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), sqrt(0.05)), length = N)^2
  j_nums <- seq(from = 0, to = R, by = 1)
  jump_mid_points <- 0
  
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
  
}

# Model from Heston (1993):
gridMaker.Heston <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  sigma <- unlist(sigma_x_params[1])
  theta <- unlist(mu_x_params[2]); kappa <- unlist(mu_x_params[1])
  var_mid_points <- seq(from = sqrt(max(theta - (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa), 0.0000001)), to = max(sqrt(theta + (3 + log(N)) * sqrt(0.5 * theta * sigma^2 / kappa)), sqrt(0.05)), length = N)^2
  j_nums <- 0
  jump_mid_points <- 0
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}

# Model from Pitt, Malik, and Doucet (2014)
gridMaker.PittMalikDoucet <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  theta <- unlist(mu_x_params[2]); phi <- unlist(mu_x_params[1])
  sigma <- unlist(sigma_x_params[1])
  
  mean <- theta
  sd <- sqrt((sigma^2) / (1 - phi^2))
  # One node point at the mean and floor(N/2) points on each side
  var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
  # Pick the first N points (or else some grids generated this way have N+1 points)
  var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
  j_nums <- seq(from = 0, to = R, by = 1)
  jump_mid_points <- 0
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}

# Model of  Taylor (1986) with leverage effect
gridMaker.TaylorWithLeverage <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  theta <- unlist(mu_x_params[2]); phi <- unlist(mu_x_params[1])
  sigma <- unlist(sigma_x_params[1])
  
  mean <- theta
  sd <- sqrt((sigma^2) / (1 - phi^2))
  
  # One node point at the mean and floor(N/2) points on each side
  var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
  # Pick the first N points (or else some grids generated this way have N+1 points)
  var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
  j_nums <- 0
  jump_mid_points <- 0
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}

# Model of  Taylor (1986)
gridMaker.Taylor <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  theta <- unlist(mu_x_params[2]); phi <- unlist(mu_x_params[1])
  sigma <- unlist(sigma_x_params[1])
  
  mean <- theta
  sd <- sqrt((sigma^2) / (1 - phi^2))
  
  # One node point at the mean and floor(N/2) points on each side
  var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
  # Pick the first N points (or else some grids generated this way have N+1 points)
  var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
  j_nums <- 0
  jump_mid_points <- 0
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}

# CAPM with SV
gridMaker.CAPM_SV <- function(model_dynamics, N, K, R, ...){ 
  mu_x_params <- model_dynamics$mu_x_params
  sigma_x_params <- model_dynamics$sigma_x_params
  theta <- unlist(mu_x_params[2]); phi <- unlist(mu_x_params[1])
  sigma <- unlist(sigma_x_params[1])
  
  mean <- theta
  sd <- sqrt((sigma^2) / (1 - phi^2))
  
  # One node point at the mean and floor(N/2) points on each side
  var_mid_points <- seq(from = 0, to = sqrt((3 + log(N)) * sd), length = floor(N / 2 + 1))^2
  # Pick the first N points (or else some grids generated this way have N+1 points)
  var_mid_points <- (sort(c(-var_mid_points[2:N], var_mid_points)) + mean)[1:N]
  j_nums <- 0
  jump_mid_points <- 0
  return(list(var_mid_points = var_mid_points, j_nums = j_nums, jump_mid_points = jump_mid_points))
}
