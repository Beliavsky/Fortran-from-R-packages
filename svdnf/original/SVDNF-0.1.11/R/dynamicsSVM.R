# a constructor function for the "SVM_dynamics" class
dynamicsSVM = function(mu = 0.038, kappa = 3.689, theta = 0.032,
                        sigma = 0.446, rho = -0.745, omega = 5.125, delta = 0.03, alpha = -0.014,
                        rho_z = -1.809, nu = 0.004, p = 0.01, phi = 0.965,
                        h = 1 / 252, coefs = NULL, model = "Heston", mu_x, mu_y, sigma_x, sigma_y,
                        jump_dist = rpois, jump_density = dpois, jump_params = 0,
                        mu_x_params, mu_y_params, sigma_x_params, sigma_y_params){
  if(kappa < 0){ 
    stop("kappa must be greater than 0")
  }
  if(sigma < 0){ 
    stop("sigma must be greater than 0")
  }
  if(rho > 1 || rho < -1) {
    stop("rho must be between -1 and 1")
  }
  if(omega < 0){
    stop("omega must be greater than 0")
  }
  if(delta < 0){
    stop("delta must be greater than 0")
  }
  if(nu < 0){ 
    stop("nu must be greater than 0")
  }
  if(p < 0 || p > 1){
    stop("p must be between 0 and 1")
  }
  if(phi < -1 || phi > 1){
    stop("phi must be between -1 and 1")
  }
  
  # Defining drifts/diffusions and jump distribution for various models (ends at line 239):
  # Model from Duffie, Pan, and Singleton (2000):
  if (model == "DuffiePanSingleton") {
    if(theta < 0){ 
      stop("theta must be greater than 0 for this model")
    }
    
    # Returns drift and diffusion
    mu_y_dps <- function(x, mu, alpha, delta, rho_z, nu, omega) {
      # Jump compensator
      alpha_bar <- exp(alpha + 0.5 * delta^2) / (1 - rho_z * nu) - 1
      
      return(h * (mu - x / 2 - alpha_bar * omega))
    }
    mu_y_params <- list(mu, alpha, delta, rho_z, nu, omega)
    sigma_y_dps <- function(x, dummy) {
      return(sqrt(h * pmax(x, 0)))
    }
    sigma_y_params <- list(0)
    
    # Volatility drift and diffusion
    mu_x_dps <- function(x, kappa, theta) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(kappa, theta)
    
    sigma_x_dps <- function(x, sigma) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma)
    
    # Jump distribution for the DuffiePanSingleton Model
    jump_density <- dpois
    jump_dist <- rpois
    jump_params <- c(h * omega)
    
    model_dynamics <- list(model = model, mu_y = mu_y_dps, sigma_y = sigma_y_dps, h = h,
      rho = rho, rho_z = rho_z, nu = nu, alpha = alpha, delta = delta,
      mu_x = mu_x_dps, sigma_x = sigma_x_dps, 
      mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
      mu_x_params = mu_x_params, sigma_x_params = sigma_x_params,
      jump_dist = jump_dist, jump_params = jump_params, jump_density = jump_density)
  }
  
  # Model from Bates (1996):
  else if (model == "Bates") {
    if(theta < 0){ 
      stop("theta must be greater than 0 for this model")
    }
    
    # Returns drift and diffusion
    mu_y_bates <- function(x, mu, alpha, delta, omega) {
      # Jump compensator
      alpha_bar <- exp(alpha + 0.5 * delta^2) - 1
      return(h * (mu - x / 2 - alpha_bar * omega))
    }
    mu_y_params <- list(mu, alpha, delta, omega)
    
    sigma_y_bates <- function(x, dummy) {
      return(sqrt(h * pmax(x, 0)))
    }
    sigma_y_params <- list(0)
    
    # Volatility drift and diffusion
    mu_x_bates <- function(x, kappa, theta) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(kappa, theta)
    
    sigma_x_bates <- function(x, sigma) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma)
    
    # Jumps for the Bates Model
    jump_density <- dpois
    jump_dist <- rpois
    jump_params <- c(h * omega)
    
    model_dynamics <- list(model = model, mu_y = mu_y_bates, sigma_y = sigma_y_bates, h = h,
      rho = rho, rho_z = 0, nu = 0, alpha = alpha, delta = delta,
      mu_x = mu_x_bates, sigma_x = sigma_x_bates,
      mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
      mu_x_params = mu_x_params, sigma_x_params = sigma_x_params,
      jump_dist = jump_dist, jump_params = jump_params, jump_density = jump_density)
  }
  
  # Model from Heston (1993):
  else if (model == "Heston") {
    if(theta < 0){ 
      stop("theta must be greater than 0 for this model")
    }
    # Returns drift and diffusion
    mu_y_heston <- function(x, mu) {
      return(h * (mu - x / 2))
    }
    mu_y_params <- list(mu)
    sigma_y_heston <- function(x, dummy) {
      return(sqrt(h * pmax(x, 0)))
    }
    sigma_y_params <- list(0)
    
    # Volatility drift and diffusion
    mu_x_heston <- function(x, kappa, theta) {
      return(x + h * kappa * (theta - pmax(0, x)))
    }
    mu_x_params <- list(kappa, theta)
    
    sigma_x_heston <- function(x, sigma) {
      return(sigma * sqrt(h * pmax(x, 0)))
    }
    sigma_x_params <- list(sigma)
    
    model_dynamics <- list(model = model, mu_y = mu_y_heston, sigma_y = sigma_y_heston, h = h,
      rho = rho, rho_z = 0, nu = 0, alpha = 0, delta = 0,
      mu_x = mu_x_heston, sigma_x = sigma_x_heston,
      mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
      mu_x_params = mu_x_params, sigma_x_params = sigma_x_params)
  }
  
  # Model from Pitt, Malik, and Doucet (2014)
  else if (model == "PittMalikDoucet") {
    
    # Returns drift and diffusion
    mu_y_pmd <- function(x, dummy) { # Dummy parameter to pass to the do.call function
      return(0)
    }
    mu_y_params <- list(0)
    sigma_y_pmd <- function(x, dummy) { # Dummy parameter to pass to the do.call function
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)
    
    # Volatility drift and diffusion
    mu_x_pmd <- function(x, phi, theta) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(phi, theta)
    sigma_x_pmd <- function(x, sigma) {
      return(sigma)
    }
    sigma_x_params <- list(sigma)
    
    # Jumps for the PMD Model
    jump_density <- dbinom
    jump_dist <- rbinom
    jump_params <- c(1, p) # size = 1, p=p for dbinom
    
    model_dynamics <- list(model = model, mu_y = mu_y_pmd, sigma_y = sigma_y_pmd,
      rho = rho, rho_z = 0, nu = 0, alpha = alpha, delta = delta,
      mu_x = mu_x_pmd, sigma_x = sigma_x_pmd,
      mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
      mu_x_params = mu_x_params, sigma_x_params = sigma_x_params,
      jump_dist = jump_dist, jump_params = jump_params, jump_density = jump_density)
  }
  
  # Model of  Taylor (1986) with the leverage effect
  else if (model == "TaylorWithLeverage") {
    
    # Returns drift and diffusion
    mu_y_twl <- function(x, dummy) { # dummy is included to pass to do.call
      return(0)
    }
    mu_y_params <- list(0)
    sigma_y_twl <- function(x, dummy) { # dummy is included to pass to do.call
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)
    
    # Volatility drift and diffusion
    mu_x_twl <- function(x, phi, theta) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(phi, theta)
    sigma_x_twl <- function(x, sigma) {
      return(sigma)
    }
    sigma_x_params <- list(sigma)
    
    model_dynamics <- list(model = model, mu_y = mu_y_twl, sigma_y = sigma_y_twl,
      rho = rho, rho_z = 0, nu = 0, alpha = 0, delta = 0,
      mu_x = mu_x_twl, sigma_x = sigma_x_twl,
      mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
      mu_x_params = mu_x_params, sigma_x_params = sigma_x_params)
  }
  # Model of  Taylor (1986)
  else if (model == "Taylor") {
    # Returns drift and diffusion
    mu_y_taylor <- function(x, dummy) { # dummy is included to pass to do.call
      return(0) # return a vector of zeroes with the length of x
    }
    mu_y_params <- list(0)
    
    sigma_y_taylor <- function(x, dummy) { # dummy is included to pass to do.call
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)
    # Volatility drift and diffusion
    mu_x_taylor <- function(x, phi, theta) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(phi, theta)
    sigma_x_taylor <- function(x, sigma) {
      return(sigma) # return a vector of sigmas with the length of x
    }
    sigma_x_params <- list(sigma)
    
    model_dynamics <- list(model = model, mu_y = mu_y_taylor, sigma_y = sigma_y_taylor,
                           rho = 0, rho_z = 0, nu = 0, alpha = 0, delta = 0,
                           mu_x = mu_x_taylor, sigma_x = sigma_x_taylor,
                           mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
                           mu_x_params = mu_x_params, sigma_x_params = sigma_x_params)
  }
  
  # Capital Asset Pricing model with stochastic volatility
  else if (model == "CAPM_SV") {
    # Default regression coefficients:
    if(is.null(coefs)){
      coefs <- c(0,1)
    }
    # Returns drift and diffusion
    mu_y_capm <- function(x, dummy) { # dummy is included to pass to do.call
      return(0) # return a vector of zeroes with the length of x
    }
    mu_y_params <- list(0)
    
    sigma_y_capm <- function(x, dummy) { # dummy is included to pass to do.call
      return(exp(x / 2))
    }
    sigma_y_params <- list(0)
    # Volatility drift and diffusion
    mu_x_capm <- function(x, phi, theta) {
      return(theta + phi * (x - theta))
    }
    mu_x_params <- list(phi, theta)
    sigma_x_capm <- function(x, sigma) {
      return(sigma) # return a vector of sigmas with the length of x
    }
    sigma_x_params <- list(sigma)
    
    model_dynamics <- list(model = model, mu_y = mu_y_capm, sigma_y = sigma_y_capm,
                           rho = 0, rho_z = 0, nu = 0, alpha = 0, delta = 0, coefs = coefs,
                           mu_x = mu_x_capm, sigma_x = sigma_x_capm,
                           mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
                           mu_x_params = mu_x_params, sigma_x_params = sigma_x_params)
  }
  else if(model == 'Custom'){
    model_dynamics <- list(model = model, mu_y = mu_y, sigma_y = sigma_y,
      rho = rho, rho_z = rho_z, nu = nu, alpha = alpha, delta = delta,
      mu_x = mu_x, sigma_x = sigma_x, coefs = coefs,
      mu_y_params = mu_y_params, sigma_y_params = sigma_y_params,
      mu_x_params = mu_x_params, sigma_x_params = sigma_x_params,
      jump_dist = jump_dist, jump_params = jump_params, jump_density = jump_density)
  }
  else{
    stop("Invalid model argument. Must select a built-in model or create custom model dynamics, see help(dynamicsSVM) for details")
  }

  class(model_dynamics) <- list('dynamicsSVM', model)
  return(model_dynamics)
}

# Print method for dynamicsSVM objects
print.dynamicsSVM <- function(x, ...){
  # Extract arguments from drift and diffusion functions
  mu_y_args <- formalArgs(x$mu_y)
  sigma_y_args <- formalArgs(x$sigma_y)
  mu_x_args <- formalArgs(x$mu_x)
  sigma_x_args <- formalArgs(x$sigma_x)
  # Remove volatility argument x and keep only model parameters
  mu_y_args <- mu_y_args[-which(mu_y_args == "x")]
  sigma_y_args <- sigma_y_args[-which(sigma_y_args == "x")]
  mu_x_args <- mu_x_args[-which(mu_x_args == "x")]
  sigma_x_args <- sigma_x_args[-which(sigma_x_args == "x")]

  if(!is.null(x$coefs)){
    k <- length(x$coefs)  
    
    coef_names <- character(k)
    
    for (i in 0:(k - 1)) {
      coef_names[i + 1] <- paste("c_", i, ":", sep = "")
    }
    
    cat("\n Regression Coefficients:\n")
    cat(rbind(coef_names, x$coefs), "\n")
  }
  
  cat("\nReturn Drift Function:\n")
  cat(deparse(x$mu_y), "\n")
  cat("\nReturn Drift Parameters:\n")
  cat(rbind(mu_y_args, unlist(x$mu_y_params)), "\n")

  
  cat("\nReturn Diffusion Function:\n")
  cat(deparse(x$sigma_y), "\n")
  cat("\nReturn Diffusion Parameters:\n")
  cat(rbind(sigma_y_args, unlist(x$sigma_y_params)), "\n")
  
  cat("\nVolatility Factor Drift Function:\n")
  cat(deparse(x$mu_x), "\n")
  cat("\nVolatility Factor Drift Parameters:\n")
  cat(rbind(mu_x_args, unlist(x$mu_x_params)), "\n")
  
  cat("\nVolatility Factor Diffusion Function:\n")
  cat(deparse(x$sigma_x), "\n")
  cat("\nVolatility Factor diffusion Parameters:\n")
  cat(rbind(sigma_x_args, unlist(x$sigma_x_params)), "\n" )
  
  if(!is.null(x$jump_density)){
  cat("\nJump Distribution\n")
  cat(deparse(x$jump_density), '\n')
  cat("\nJump Distribution Parameters\n")
  cat(rbind(formalArgs(x$jump_density)[c(2,(length(x$jump_params)+1))],
            unlist(x$jump_params)), "\n" )
  }
}

