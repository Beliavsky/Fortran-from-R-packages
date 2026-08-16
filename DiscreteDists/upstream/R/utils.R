#' logLik function for hyper Poisson
#' @description Calculates logLik for hyper Poisson distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_HYPERPO <- function(logparam=c(0, 0), x){
  return(sum(dHYPERPO(x     = x,
                      mu    = exp(logparam[1]),
                      sigma = exp(logparam[2]),
                      log=TRUE)))
}
#' Initial values for hyper Poisson
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_sigma_HYPERPO <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_HYPERPO,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' Auxiliar function for F11
#' @description This function is used inside F11 function.
#' @param x vector
#' @param tol this is the tolerance of the infinite sum.
#' @return returns a logical value if the tolerance level is met.
#' @keywords internal
#' @keywords noRd
#' @export
stopping <- function (x, tol) {
  all(abs(x) <= tol, na.rm = TRUE)
}
#' Auxiliar function for hyper Poisson
#' @description This function is used inside density function of Hyper Poisson.
#' @param z,c values for F11.
#' @param maxiter_series maximum value to obtain F11.
#' @param tol this is the tolerance of the infinite sum.
#' @return returns the value for the F11 function.
#' @keywords internal
#' @keywords noRd
#' @export
F11 <- function(z, c, maxiter_series = 10000, tol = 1.0e-10) {
  fac  <- 1
  temp <- 1
  L    <- c
  for (n in seq_len(maxiter_series)) {
    fac    <- fac * z / L
    series <- temp + fac
    if (stopping(series - temp, tol)){
      return(Re(series))
    }
    temp   <- series
    L      <- L + 1
  }
  if (tol >= 0)
    return(Re(series))
}
F11 <- Vectorize(F11)
#' Auxiliar function for hyper Poisson
#' @description This function is used to calculate (a)r.
#' @param a first value.
#' @param r second value.
#' @return returns the value for the a(r) function.
#' @keywords internal
#' @keywords noRd
#' @export
AR <- function(a, r) {
  res <- gamma(a+r) / gamma(a)
  res
}
#' Auxiliar function to generate values for hyper Poisson
#' @description This function is used inside random function of Hyper Poisson.
#' @param sigma value for sigma parameter.
#' @param mu value for mu parameter.
#' @keywords internal
#' @keywords noRd
#' @export
simulate_hp <- function(sigma, mu) {
  pochammer <- function(a, r) if (r == 0) 1 else prod(a:(a + r - 1))
  u <- runif(1)
  y <- 0
  p <- 0
  value <- F11(mu, sigma)
  while (p < u) {
    p <- p + mu ^ y / (value * pochammer(sigma, y))
    y <- y + 1
  }
  y - 1
}
#' logLik function for hyper Poisson in second parameterization
#' @description Calculates logLik for hyper Poisson distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_HYPERPO2 <- function(logparam=c(0, 0), x){
  return(sum(dHYPERPO2(x     = x,
                       mu    = exp(logparam[1]),
                       sigma = exp(logparam[2]),
                       log=TRUE)))
}
#' Initial values for hyper Poisson in second parameterization
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_sigma_HYPERPO2 <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_HYPERPO2,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' Auxiliar function to obtain lambda from E(X) in HYPERPO2
#' @description This function implements the procedure given in page 152.
#' @param media the value for the mean or E(X).
#' @param gamma the value for the gamma parameter.
#' @return returns the value of lambda to ensure the mean and gamma.
#' @keywords internal
#' @keywords noRd
#' @export
obtaining_lambda <- function(media, gamma) {
  # Begin aux function, based on expression 6 of Saez-Castillo (2013)
  fun_exp_6 <- function(x, media, gamma) x-(gamma-1)*(1-1/f11_cpp(gamma, x))-media
  fun_exp_6 <- Vectorize(fun_exp_6)
  # End aux function
  if (gamma == 1)
    result <- media
  else {
    mini <- min(media, max(media+gamma-1, gamma*media))
    maxi <- max(media, min(media+gamma-1, gamma*media))
    res <- try(uniroot(f=fun_exp_6_vec_cpp,
                       #f=fun_exp_6,
                       lower=mini, upper=maxi,
                       media=media, gamma=gamma))
    result <- ifelse(class(res)=="try-error", media, res$root)
  }
  result
}
obtaining_lambda <- Vectorize(obtaining_lambda)
#' logLik function for Discrete Burr Hatke
#' @description Calculates logLik for Discrete Burr Hatke  distribution.
#' @param param value for mu.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_DBH <- function(param=0.5, x){
  mu <- param
  minus_ll <- -sum(dDBH(x=x, mu=mu, log=TRUE))
  return(minus_ll)
}
#' Initial values for Discrete Burr Hatke
#' @description This function generates initial values for the parameter mu.
#' @param y vector with the response variable.
#' @return returns a scalar with the MLE estimation.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optimize
estim_mu_DBH <- function(y){
  mod <- optimize(f=logLik_DBH, interval=c(0, 1), x=y)
  res <- mod$minimum
  return(res)
}
#' logLik function for Discrete Lindley distribution
#' @description Calculates logLik for Discrete Lindley  distribution.
#' @param param value for mu.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_DLD <- function(param=0.5, x){
  mu <- param
  minus_ll <- -sum(dDLD(x=x, mu=mu, log=TRUE))
  return(minus_ll)
}
#' Initial values for Discrete Lindley
#' @description This function generates initial values for the parameter mu.
#' @param y vector with the response variable.
#' @return returns a scalar with the MLE estimation.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_DLD <- function(y){
  res1 <- optim(par=0.5, fn=logLik_DLD, method='L-BFGS-B',
                lower=1e-10, upper=Inf, x = y)
  res <- res1$par
  return(res)
}

#' logLik function for discrete Inverted Kumaraswamy
#' @description Calculates logLik for discrete Inverted Kumaraswamy distribution.
#' @param param vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_DIKUM <- function(param=c(0, 0), x){
  return(sum(dDIKUM(x = x,
                   mu = exp(param[1]),
                   sigma = exp(param[2]),
                   log=TRUE)))
}
#' Initial values for discrete Inverted Kumaraswamy
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_sigma_DIKUM <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_DIKUM,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' logLik function for Poisson XLindley distribution
#' @description Calculates logLik for Poisson XLindley distribution distribution.
#' @param param parameter mu in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_POISXL <- function(param=0, x){
  return(sum(dPOISXL(x = x,
                     mu = exp(param),
                     log=TRUE)))
}
#' Initial values for discrete Poisson XLindley distribution
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a scalar with the MLE estimation.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_POISXL <- function(y) {
  mod <- optim(par=c(0),
               fn=logLik_POISXL,
               method="Brent",
               control=list(fnscale=-1, maxit=100000),
               x=y,
               lower=-100, upper=100)
  res <- exp(mod$par)
  names(res) <- c("mu_hat")
  return(res)
}
#' logLik function for GGEO
#' @description Calculates logLik for GGEO distribution.
#' @param param vector with parameters in log and logit scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_GGEO <- function(param=c(0, 0), x){
  inv_logit <- function(x) 1/(1 + exp(-x))
  return(sum(dGGEO(x,
                   mu    = inv_logit(param[1]),
                   sigma = exp(param[2]),
                   log=TRUE)))
}
#' Initial values for GGEO
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_sigma_GGEO <- function(y) {
  inv_logit <- function(x) 1/(1 + exp(-x))
  mod <- optim(par=c(0, 0),
               fn=logLik_GGEO,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = inv_logit(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' logLik function for DGEII
#' @description Calculates logLik for DGEII distribution.
#' @param transf_param vector with parameters in log and logit scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @keywords noRd
#' @export
logLik_DGEII <- function(transf_param=c(0, 0), x){
  inv_logit <- function(x) 1/(1 + exp(-x))
  return(sum(dDGEII(x,
                    mu    = inv_logit(transf_param[1]),
                    sigma = exp(transf_param[2]),
                    log=TRUE)))
}
#' Initial values for DGEII
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_sigma_DGEII <- function(y) {
  inv_logit <- function(x) 1/(1 + exp(-x))
  mod <- optim(par=c(1 - 1/(1+mean(y)), 0),
               fn=logLik_DGEII,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = inv_logit(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
