# Run help(pars) for description and argument details.
pars <- function(dynamics, ...) UseMethod("pars") 

pars.default <- function(dynamics, ...){
  stop("This class of state-space models is not supported by the pars function.")
}

pars.dynamicsSVM <- function(dynamics, rho = NULL, delta = NULL, alpha = NULL, rho_z = NULL, nu = NULL, jump_params_list = "dummy", ...){
  # vector of commom model parameters
  par_v <- c(rho, delta, alpha, rho_z, nu)
  # Check input parameters
  if(!(is.null(par_v) | all(par_v =='var'))){
    stop(paste0("Parameters (rho, delta, alpha, and rho_z) should be set to NULL if they is not in the model. Otherwise, it should be set to 'var'."))
  }
  
  model <- dynamics$model
  
  coefs <- dynamics$coefs
  
  if (model == "DuffiePanSingleton") {
    rho <- 'var'
    alpha <- 'var'
    delta <- 'var'
    nu <- 'var'
    rho_z <- 'var'
  }
  # Model from Bates (1996):
  else if (model == "Bates") {
    rho <- 'var'
    alpha <- 'var'
    delta <- 'var'
  }
  
  # Model from Heston (1993):
  else if (model == "Heston") {
    rho <- 'var'
  }
  
  # Model from Pitt, Malik, and Doucet (2014)
  else if (model == "PittMalikDoucet") {
    rho <- 'var'
    alpha <- 'var'
    delta <- 'var'
    
    jump_params_list <- "p"
    
  }
  
  # Model of  Taylor (1986) with leverage effect
  else if  (model == "TaylorWithLeverage") {
    rho <- 'var'
  }

  # Making a list of the parameter names in the optimization: 
  
  mu_y <- dynamics$mu_y; sigma_y <- dynamics$sigma_y
  mu_x <- dynamics$mu_x; sigma_x <- dynamics$sigma_x
  
  # Getting the arguments from drift/diffusion functions
  args_vec <- c(formalArgs(mu_y), formalArgs(sigma_y),
                formalArgs(mu_x), formalArgs(sigma_x))
  
  # Adding possible parameters that aren't in the drift/diffusion
  # If set to "var", they will be optimized, or else they are fixed at given values.
  if(!is.null(rho)){
  if(rho == 'var'){
    args_vec <- c(args_vec, 'rho')
  }
  }
  if(!is.null(delta)){
  if(delta == 'var'){
    args_vec <- c(args_vec, 'delta')
  } 
  }
  if(!is.null(alpha)){
  if(alpha == 'var'){
    args_vec <- c(args_vec, 'alpha')
  }
  }
  if(!is.null(rho_z)){
  if(rho_z == 'var'){
    args_vec <- c(args_vec, 'rho_z')
  }
  }
  if(!is.null(nu)){
  if(nu == 'var'){
    args_vec <- c(args_vec, 'nu')
  }
}
  
  if(!is.null(coefs)){
    k <- length(coefs)  
    
    coef_names <- character(k)
    
    for (i in 0:(k - 1)) {
      coef_names[i + 1] <- paste("c_", i, sep = "")
    }
    args_vec <- c(coef_names, args_vec)
  }
  args_vec <- c(args_vec, jump_params_list)
  
  args_vec <- args_vec[-which(args_vec == "dummy")]
  args_vec <- args_vec[-which(args_vec == "x")]
  args_vec <- unique(args_vec)
  return(args_vec)
}
