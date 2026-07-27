modelSim <- function(dynamics, ...) UseMethod("modelSim") 

modelSim.default <- function(dynamics, ...){
  stop("This class of state-space models is not supported by the DNF function.")
}
# Run help(modelSim) for description and argument details.
modelSim.dynamicsSVM <- function(dynamics, t, init_vol = 0.032, ...) {
  # Default jump vectors
  j_x_vec <- rep(0, length = t)
  j_y_vec <- rep(0, length = t)
  j_nums <- rep(0, length = t)
  
  if (! is.null(dynamics$jump_dist)) { # Generating jumps if there is a jump distribution in the model
    j_nums <- do.call(dynamics$jump_dist, c(list(t), dynamics$jump_params))
    j_x_vec <- rgamma(t, shape = j_nums, scale = dynamics$nu)
    j_y_vec <- rnorm(t, mean = j_nums *  dynamics$alpha, sd = sqrt(j_nums * dynamics$delta^2)) + dynamics$rho_z * j_x_vec
  }
  
  mu_y <- dynamics$mu_y; sigma_y <- dynamics$sigma_y
  mu_y_params <- dynamics$mu_y_params; sigma_y_params <- dynamics$sigma_y_params
  
  mu_x <- dynamics$mu_x; sigma_x <- dynamics$sigma_x
  mu_x_params <- dynamics$mu_x_params; sigma_x_params <- dynamics$sigma_x_params
  
  rho = dynamics$rho
  # Draw our noise
  epsilon_x <- rnorm(t, 0, 1)
  epsilon_y <- rnorm(t, 0, 1)
  x_zero <- init_vol # Assume x_0 = theta
  x_vec <- c(x_zero)
  y_vec <- c()
  for (i in 1:t) {
    x_tmin1 <- x_vec[i]
    
    # Evaluate the drift/diffusion:
    mu_y_eval <- do.call(mu_y, c(list(x_tmin1), mu_y_params))
    mu_x_eval <- do.call(mu_x, c(list(x_tmin1), mu_x_params))
    
    sigma_y_eval <- do.call(sigma_y, c(list(x_tmin1), sigma_y_params))
    sigma_x_eval <- do.call(sigma_x, c(list(x_tmin1), sigma_x_params))
    
    # Generate returns and the volatility factor:
    y_t <- mu_y_eval + sigma_y_eval * (rho * epsilon_x[i] + sqrt(1 - rho^2) * epsilon_y[i]) + j_y_vec[i]
    x_t <- mu_x_eval + sigma_x_eval * epsilon_x[i] + j_x_vec[i]
    
    # Update the sequences:
    x_vec <- c(x_vec, x_t)
    y_vec <- c(y_vec, y_t)
  }
  x_vec <- x_vec[-1] # Remove x_zero = theta
  model_sim <- list(volatility_factor = x_vec, returns = y_vec)
  return(model_sim)
}

