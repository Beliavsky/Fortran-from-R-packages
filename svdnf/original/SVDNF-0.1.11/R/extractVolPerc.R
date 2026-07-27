# Function to find the percentile p of a filtering/prediction distribution grid.
# Run help(DNF) for description and argument details.
extractVolPerc <- function(x, ...) UseMethod("extractVolPerc") 

extractVolPerc.default <- function(x, ...){
  stop("This class of objects is not supported by the extractVolPerc function.")
}
extractVolPerc.SVDNF<- function(x, p = 0.5, pred = F, ...) {
  # Extract grids from x
  var_mid_points <- x$grids$var_mid_points
  jump_mid_points <- x$grids$jump_mid_points
  j_num <- x$grids$j_nums
  
  # Extract dynamics from x
  dynamics = x$dynamics
  # Define N, K, and R from grids
  N <- length(var_mid_points); K <- length(jump_mid_points); R <- max(j_num)
  
  # Filtering distribution from the DNF and length of the series
  filtering <- x$filter_grid
  T <- ncol(filtering)-1 # remove one as filtering from t = 0 ,..., t = T
  
  NNKR_grid <- expand.grid(var_mid_points,var_mid_points, jump_mid_points, j_num)
  v_t = unlist(NNKR_grid[1]); v_tmin1 <- unlist(NNKR_grid[2]); j_mt <- unlist(NNKR_grid[3]); j_nums <- unlist(NNKR_grid[4])
  
  # Transition probabilities : 
  p_v <- dnorm(v_t, do.call(dynamics$mu_x, c(list(v_tmin1), dynamics$mu_x_params)), sd = do.call(dynamics$sigma_x, c(list(v_tmin1), dynamics$sigma_x_params)))
  
  p_n <- 1 # default value of 1, if there are no jumps in the model
  p_j_vol <- 1 # default value of 1, if there are no volatility jumps in the model
  
  if (any(j_nums != 0)) { # If there are no jumps, leave jump_mat = 1.
    # Probability of having j_nums jumps
    p_n <- do.call(dynamics$jump_density, c(list(j_nums), dynamics$jump_params))
  }
  
  if (any(jump_mid_points != 0)) {
    # Probability of having jumps of certain size within the jump interval grid.
    p_j_vol <- dgamma(j_mt, shape = j_nums, scale = dynamics$nu)
  }
  
  Transition <- p_v * p_n * p_j_vol
  #Sum across jump size:
  Transition <- rowSums(matrix(Transition, nrow = N * N * (R + 1), ncol = K, byrow = F))
  #Sum across jump number:
  Transition <- rowSums(matrix(Transition, nrow = N * N, ncol = R + 1, byrow = F))
  Transition <- matrix(Transition, nrow = N, ncol = N) 
  
  pred_mat <- matrix(NA, nrow = N, ncol = T)
  # Compute the prediction distribution at each step 
  for(i in (1:T)){
    pred_mat[N:1,i] <- (x$filter_grid[N:1,i]) %*% t(Transition)
  }
  
  filtering_CDF <- filtering
  pred_CDF <- pred_mat
  
  # Compute the volatility's CDF at each time point from the filtering distribution
  for (i in (2:N)) {
    filtering_CDF[i, ] <- filtering_CDF[(i - 1), ] + filtering[i, ]
    pred_CDF[i, ] <- pred_CDF[(i - 1), ] + pred_mat[i, ]
  }
  # Normalize so the CDF sums to 1
  filtering_CDF <- filtering_CDF %*% diag(1 / filtering_CDF[N, ])
  pred_CDF <- pred_CDF %*% diag(1 / pred_CDF[N, ])
  
  CDF <- filtering_CDF
  if(pred){
    CDF <- pred_CDF
  }
  per_vec <- c()
  T <- ncol(CDF)
  N <- nrow(CDF)
  p <- 1 - p
  for (i in (1:T)) {
    per_vec <- c(per_vec, var_mid_points[N - which((CDF[, i]) > p)[1] + 1])
  }
  return(per_vec)
}

extractVolPerc.DNFOptim<- function(x, p = 0.5, pred = F, ...) {
  dnf <- x$SVDNF
  per_vec <- extractVolPerc(dnf, p, pred, ...)
  return(per_vec)
}
