# Function probCalculator:
#   Computes and multiplies together all the probabilities in Equation (2)
#   of the SVDNF package vignette besides the t-1 filtering distribution.

# Inputs:
#    - grids : List of nodes to use for numerical integration.
#    - mu_x, mu_y, sigma_x, sigma_y : Drift and diffusion functions
#                                     used for the DNF (see help(DNF) for a description).
#    - mu_x_params, mu_y_params, sigma_x_params, sigma_y_params : Lists of parameters to be
#                                                                 passed in the drift/diffusion
#                                                                 function
#    - model : Model for which we are applying the DNF.
# Output:
#   	- d_probs : A N*N*K*(R+1) x t matrix containing the product of the probabilities in Equation (2)
#                 to be multiplied by the previous time's filtering distribution.

probCalculator <- function(grids, R = 1, N = 50, K = 20, data, dynamics) {
  var_mid_points <- grids$var_mid_points
  j_nums <- grids$j_nums
  jump_mid_points <- grids$jump_mid_points

  # Interval vectors
  var_intervals <- c(var_mid_points[1] - (var_mid_points[2] - var_mid_points[1]), var_mid_points, Inf)
  var_intervals <- (var_intervals[1:(N + 1)] + var_intervals[2:(N + 2)]) / 2

  jump_intervals <- c(-jump_mid_points[1], jump_mid_points, Inf)
  jump_intervals <- (jump_intervals[1:(K + 1)] + jump_intervals[2:(K + 2)]) / 2

  # Using expand.grid to get all combination of our grids
  NNKR_grid <- expand.grid(var_mid_points, var_mid_points, jump_mid_points, j_nums)

  x_t <- unlist(NNKR_grid[, 1])
  x_tmin1 <- unlist(NNKR_grid[, 2])
  j_mt <- unlist(NNKR_grid[, 3])
  j_nums <- unlist(NNKR_grid[, 4])

  # Evaluating the drift/diffusion at the mid points
  mu_y_eval <- do.call(dynamics$mu_y, c(list(x_tmin1), dynamics$mu_y_params))
  mu_x_eval <- do.call(dynamics$mu_x, c(list(x_tmin1), dynamics$mu_x_params))

  sigma_y_eval <- do.call(dynamics$sigma_y, c(list(x_tmin1), dynamics$sigma_y_params))
  sigma_x_eval <- do.call(dynamics$sigma_x, c(list(x_tmin1), dynamics$sigma_x_params))

  # Mean and variance of y conditional on the latent factors
  eps <- (x_t - mu_x_eval - j_mt) / (sigma_x_eval)

  mus <- mu_y_eval + dynamics$rho * sigma_y_eval * eps + dynamics$alpha * j_nums + dynamics$rho_z * j_mt
  sigmas <- sqrt((1 - dynamics$rho^2) * (sigma_y_eval^2) + j_nums * dynamics$delta^2)
  # Mean and variance of x_t given x_t_min1, j_t^x, and n_t
  mu_v <- mu_x_eval + j_mt
  sig_v <- sigma_x_eval

  # Computing the probability of x_t being within certain intervals
  q_v <- pnorm(var_intervals[findInterval(x_t, vec = var_intervals) + 1], mu_v, sig_v) - pnorm(var_intervals[findInterval(x_t, vec = var_intervals)], mu_v, sig_v)
  p_n <- 1 # default value of 1, if there are no jumps in the model
  q_j <- 1 # default value of 1, if there are no volatility jumps in the model

  if (any(j_nums != 0)) { # If there are no jumps, leave jump_mat = 1.
    # Probability of having j_nums jumps
    p_n <- do.call(dynamics$jump_density, c(list(j_nums), dynamics$jump_params))
  }

  if (any(jump_mid_points != 0)) {
    # Proability of having jumps of certain size within the jump interval grid.
    q_j <- pgamma(jump_intervals[findInterval(j_mt, vec = jump_intervals) + 1], shape = j_nums, scale = dynamics$nu) - pgamma(jump_intervals[findInterval(j_mt, vec = jump_intervals)], shape = j_nums, scale = dynamics$nu)
  }
  jump_mat <- p_n * q_j

  d_probs <- dnorm_cpp_prod(data, mus, sigmas, jump_mat * q_v) # C++ function to compute measurement equation times the other probabilities.
  return(d_probs)
}
