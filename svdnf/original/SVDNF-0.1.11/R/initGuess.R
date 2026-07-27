# Run help(initGuess) for description and argument details.
initGuess <- function(dynamics, ...) UseMethod("initGuess") 

initGuess.default <- function(dynamics, ...){
  stop("This class of state-space models is not supported by the initGuess function.")
}

initGuess.dynamicsSVM <- function(dynamics, data, factors = NULL, N = 50, K = 20, R = 1, grids = "Default", ...){
  # Extract returns if xts object.
  if(is.xts(data)){
    # Extract returns to run DNF with
    ret <- as.vector(coredata(data))
  }else{
    ret <- data
  }
  model <- dynamics$model
  # Check if a built-in model was given.
  if(!any(model %in% c("Taylor", "TaylorWithLeverage", "PittMalikDoucet", "Heston", "Bates", "DuffiePanSingleton", "CAPM_SV"))){
    stop("The model used is not supported by the initGuess function, please use a built-in model or pass initial parameters to the DNFOptim function using the par argument.")
  }
  T_max <- length(ret)
  # Fit exponentially weighted moving average
  lambda <- 0.9 
  # t = 0 var
  EWMA_vol <- var(ret)
  for(i in (1:T_max)){
    EWMA_vol <- c(EWMA_vol, EWMA_vol[i] *lambda +(1-lambda)*ret[i]^2)
  }
  
  X_t <- EWMA_vol[2:(T_max+1)]
  X_tmin1 <- EWMA_vol[1:T_max]
  
  if(any(model %in% c("Taylor", "TaylorWithLeverage", "PittMalikDoucet", "CAPM_SV"))){
    theta <- log(var(ret))
    dynamics$mu_x_params[2] <- theta
    
    vol_df <- data.frame(y = log(X_t), X = log(X_tmin1))
    lin_reg <- lm(y ~ X, vol_df)
    phi <- max(min(lin_reg$coefficients[2], 0.99), -0.99)
    #theta <- lin_reg$coefficients[1] / (1-phi)
    sigma <- max(summary(lin_reg)$sigma, 0.1)
    par <- c(phi, theta, sigma)
    dynamics <- dynamicsSVM(model = model, phi = phi, theta = theta, sigma = sigma)
    if(model == "CAPM_SV"){
      par <- c(dynamics$coefs, par)
    }
    if((model != "Taylor") & (model != "CAPM_SV")){
      # Estimate rho
      epsilon_y <- ret / (exp(X_t / 2))
      epsilon_x <- (X_t - (theta + phi * (X_tmin1 - theta))) /sigma
      rho <- cor(epsilon_x, epsilon_y)
      
      #dynamics$rho <- rho
      par <- c(par, rho)
      dynamics <- dynamicsSVM(model = model, phi = phi, theta = theta, sigma = sigma, rho = rho)
    }
    
  }
  if(model == "CAPM_SV"){
    lin_reg <- lm.fit(y = ret, x = t(factors))
    dynamics$coefs[1] <- lin_reg$coefficients[1]
    dynamics$coefs[2] <- lin_reg$coefficients[2]
  }
  
  if(any(model %in% c("Heston", "Bates", "DuffiePanSingleton"))){
    h <- dynamics$h
    theta <- var(ret)/h
    mu <- mean(ret)/h
    vol_df <- data.frame(y = X_t, X = X_tmin1 )
    lin_reg <- lm(y ~ X, vol_df)
    kappa <- (1 - lin_reg$coefficients[2]) / h
    sigma <- summary(lin_reg)$sigma / (sqrt(theta * h))
    # Estimate rho
    epsilon_y <- (ret - mean(ret)) / sqrt(X_t * h)
    epsilon_x <- (X_t - (X_tmin1 + h * kappa * (theta - X_tmin1)))/ (sigma*sqrt(X_t*h))
      
    # Update dynamics parameters
    rho <- cor(epsilon_x, epsilon_y)
    dynamics$rho <- cor(epsilon_x, epsilon_y)
      
    par <-c(mu, kappa, theta, sigma, rho)
    dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h)
    }
  # Run initial DNF
  init_DNF <- DNF(data = ret, dynamics = dynamics, factors = factors, N = N, K = K, R = R)
  # init while loop
  i <- 1
  # Store likelihoods to check for improvements
  likelihoods <- c(0, init_DNF$log_likelihood)
  filtered_mean <- colSums(sort(init_DNF$grids$var_mid_points, decreasing = T) * init_DNF$filter_grid)
  X_t <- filtered_mean[2:(T_max+1)]
  X_tmin1 <- filtered_mean[1:T_max]
  
  if(any(model %in% c("Taylor", "TaylorWithLeverage", "PittMalikDoucet", "CAPM_SV"))){
    vol_df <- data.frame(y = X_t, X = X_tmin1)
    lin_reg <- lm(y ~ X, vol_df)
    phi <- max(min(lin_reg$coefficients[2], 0.99), -0.99)
    sigma <- max(summary(lin_reg)$sigma, 0.1)
    par <- c(phi, theta, sigma)
    dynamics <- dynamicsSVM(model = model, phi = phi, theta = theta, sigma = sigma)
    if(model == "CAPM_SV"){
      par <- c(dynamics$coefs, par)
    }
    if((model != "Taylor") & (model != "CAPM_SV")){
      # Estimate rho
      epsilon_y <- ret / (exp(X_t / 2))
      epsilon_x <- (X_t - (theta + phi * (X_tmin1 - theta))) /sigma
      rho <- cor(epsilon_x, epsilon_y)
      
      par <- c(par, rho)
      dynamics <- dynamicsSVM(model = model, phi = phi, theta = theta, sigma = sigma, rho = rho)
    }
    
  }
  
  if(any(model %in% c("Heston", "Bates", "DuffiePanSingleton"))){
    vol_df <- data.frame(y = X_t, X = X_tmin1)
    lin_reg <- lm(y ~ X, vol_df)
    kappa <- (1 - lin_reg$coefficients[2]) / h
    sigma <- summary(lin_reg)$sigma / (sqrt(theta * h))
    # Estimate rho
    epsilon_y <- (ret - mean(ret)) / sqrt(X_t * h)
    epsilon_x <- (X_t - (X_tmin1 + h * kappa * (theta - X_tmin1)))/ (sigma*sqrt(X_t*h))
    
    # Update dynamics parameters
    rho <- cor(epsilon_x, epsilon_y)
    dynamics$rho <- cor(epsilon_x, epsilon_y)
    
    par <-c(mu, kappa, theta, sigma, rho)
    dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h)
  }
  
  # Jump parameters
  if(any(model %in% c("PittMalikDoucet", "Bates", "DuffiePanSingleton"))){
    if(any(model %in% c("Bates", "DuffiePanSingleton"))){
      # 99% CI to detect outlier returns
      ucl <- mu * h + qnorm(0.995) * sqrt(h * X_t)
      lcl <- mu * h - qnorm(0.995) * sqrt(h * X_t)
    } else{
      ucl <- qnorm(0.995) * exp(X_t / 2)
      lcl <- - qnorm(0.995) * exp(X_t / 2) 
    }
    jump_test_df <- data.frame(ret = ret, ucl = ucl, lcl = lcl) 
    # Which returns fall outside of the 99% CI
    outliers_ind <- (ret<= lcl | ret >= ucl)
    outliers <- jump_test_df[ret<= lcl | ret >= ucl,]
    jump_prop <- nrow(outliers)/T_max
    alpha <- mean(outliers$ret)
    delta <- sd(outliers$ret)
    if(is.na(delta) |is.na(alpha)){
      #set default
      alpha = -0.01
      delta = 0.01
    }
    # Excess large returns counted as jumps
    if(model == "PittMalikDoucet"){
      p  <- max(jump_prop - 0.01, 0)
      #dynamics$jump_params <- c(1, p)
      par <- c(par,  delta, alpha, p)
      dynamics <- dynamicsSVM(model = model,theta = theta, sigma = sigma, rho = rho, alpha = alpha, delta = delta, phi = phi, p = p)
    }else{
      omega <- max((jump_prop - 0.01) /h, 0)
      dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h, alpha = alpha, delta = delta, omega = omega)
      par <- c(mu, alpha, delta, omega, kappa, theta, sigma, rho)
      
    }
    if(model == "DuffiePanSingleton"){ 
      mean_vol <- X_tmin1 + h * kappa * (theta - X_tmin1)
      nu <- max(mean((X_t -  mean_vol)[outliers_ind]),0.005) 
      par <- c(par[1:3], 0, nu, par[4:8])
      dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h, alpha = alpha, delta = delta, omega = omega, nu = nu)
    }
  }
  par_mat <- par
    while(i < 20){
     run_DNF <- DNF(data = ret, dynamics = dynamics, factors = factors, N = N, K = K, R = R, grids = grids)
   likelihoods <- c(likelihoods, run_DNF$log_likelihood)
   
   filtered_mean <- colSums(sort(run_DNF$grids$var_mid_points, decreasing = T) * run_DNF$filter_grid)
   X_t <- filtered_mean[2:(T_max+1)]
   X_tmin1 <- filtered_mean[1:T_max]
   
   if(any(model %in% c("Taylor", "TaylorWithLeverage", "PittMalikDoucet", "CAPM_SV"))){
     vol_df <- data.frame(y = X_t, X = X_tmin1)
     lin_reg <- lm(y ~ X, vol_df)
     phi <- max(min(lin_reg$coefficients[2], 0.99), -0.99)
     sigma <- max(summary(lin_reg)$sigma, 0.1)

     
     par <- c(phi, theta, sigma)
     dynamics <- dynamicsSVM(model = model, phi = phi, theta = theta, sigma = sigma)
     if(model == "CAPM_SV"){
       par <- c(dynamics$coefs, par)
     }
     if((model != "Taylor") & (model != "CAPM_SV")){
       # Estimate rho
       epsilon_y <- ret / (exp(X_t / 2))
       epsilon_x <- (X_t-  (theta + phi * (X_tmin1 - theta))) /sigma
       rho <- cor(epsilon_x, epsilon_y)
       
       par <- c(par, rho)
       dynamics <- dynamicsSVM(model = model, phi = phi, theta = theta, sigma = sigma, rho = rho)
     }

   }
   
   if(any(model %in% c("Heston", "Bates", "DuffiePanSingleton"))){
     vol_df <- data.frame(y = X_t, X = X_tmin1)
     lin_reg <- lm(y ~ X, vol_df)
     kappa <- (1 - lin_reg$coefficients[2]) / h
     sigma <- summary(lin_reg)$sigma / (sqrt(theta * h))
     # Estimate rho
     epsilon_y <- (ret - mean(ret)) / sqrt(X_t * h)
     epsilon_x <- (X_t - (X_tmin1 + h * kappa * (theta - X_tmin1)))/ (sigma*sqrt(X_t*h))
     
     # Update dynamics parameters
     rho <- cor(epsilon_x, epsilon_y)
    
     par <-c(mu, kappa, theta, sigma, rho)
     dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h)
   }
   
   # Jump parameters
   if(any(model %in% c("PittMalikDoucet", "Bates", "DuffiePanSingleton"))){
     if(any(model %in% c("Bates", "DuffiePanSingleton"))){
       # 99% CI to detect outlier returns 
       ucl <- mu * h + qnorm(0.995) * sqrt(h * X_t)
       lcl <- mu * h - qnorm(0.995) * sqrt(h * X_t)
     } else{
       ucl <- qnorm(0.995) * exp(X_t / 2)
       lcl <- - qnorm(0.995) * exp(X_t / 2) 
     }
     jump_test_df <- data.frame(ret = ret, ucl = ucl, lcl = lcl) 
     # Which returns fall outside of the 99% CI
     outliers_ind <- (ret<= lcl | ret >= ucl)
     outliers <- jump_test_df[ret<= lcl | ret >= ucl,]
     jump_prop <- nrow(outliers) / T_max

     alpha <- mean(outliers$ret)
     delta <- sd(outliers$ret)
     if(is.na(delta) |is.na(alpha)){
       #set default
       alpha = -0.01
       delta = 0.01
     }
     # Excess large returns counted as jumps
     if(model == "PittMalikDoucet"){
       p  <- max(jump_prop - 0.01, 0)
       par <- c(par,  delta, alpha, p)
       dynamics <- dynamicsSVM(model = model,theta = theta, sigma = sigma, rho = rho, alpha = alpha, delta = delta, phi = phi, p = p)
     }else{
       omega <- max((jump_prop - 0.01) /h, 0)
       dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h, alpha = alpha, delta = delta, omega = omega)
       par <- c(mu, alpha, delta, omega, kappa, theta, sigma, rho)
       
     }
     if(model == "DuffiePanSingleton"){ 
       mean_vol <- X_tmin1+h*kappa*(theta-X_tmin1)
       nu <- max(mean((X_t -  mean_vol)[outliers_ind]),0.005)
       par <- c(par[1:3], 0, nu, par[4:8])
       dynamics <- dynamicsSVM(model = model, mu = mu, kappa = kappa, theta = theta, sigma = sigma, rho = rho, h = h, alpha = alpha, delta = delta, omega = omega, nu = nu)
     }
   }
   par_mat <- rbind(par_mat, par)
    i <- i +1 
    }
  return(par_mat[which.max(likelihoods[-1]),])
}
