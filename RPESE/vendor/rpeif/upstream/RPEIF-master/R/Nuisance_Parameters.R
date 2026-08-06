#' @title Nuisance Parameters Computation
#' 
#' @description \code{nuis.pars} returns the value of the nuisance parameters used in the evaluation of the shape of influence functions for risk and performance measures.
#'
#' @param mu Mean parameter.
#' @param sd Standard deviation parameter.
#' @param c Constant value for threshold.
#' @param alpha Parameters for the lower tail quantile.
#' @param beta Parameter for the upper tail quantile.
#'
#' @return List of nuisance parameters.
#' 
#' @details 
#' For further details on the usage of the \code{nuisParsFn} function, please refer to Section 3.1 for the \code{RPEIF} vignette.
#'
#' @author Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
#'
#' @export
#'
#' @examples
#' # Nuisance parameters using default values
#' defaultNuisance <- nuisParsFn()
#'
#' # Nuisance parameters using specified values
#' specifiedNuisance <- nuisParsFn(mu=0.02, sd=0.1, c=0.01, alpha=0.05, beta=0.1)
#'
nuisParsFn  <- function(mu=0.01, sd=0.05, c=0, alpha=0.1, beta=0.1){
  
  # Check format of input
  if(!is.numeric(c(mu, sd, c, alpha, beta)))
    stop("All function arguments must be numeric.")
  
  # Check input for alpha parameter
  if(alpha < 0 || alpha >1)
    stop("The argument alpha must be a numeric value between 0 and 1.")
  
  # Check input for beta parameter
  if(beta < 0 || beta >1)
    stop("The argument beta must be a numeric value between 0 and 1.") 
    
  # Nuisance Parameters for Risk Estimators
  
  # ES
  q.alpha <- mu+qnorm(alpha)*sd
  es.alpha <- -mu+dnorm(qnorm(alpha))/alpha*sd
  
  # LPM1
  d <- (c-mu)/sd
  lpm1 <- (d*pnorm(d)+dnorm(d))*sd
  
  # LPM2
  lpm2 <- ((d^2+1)*pnorm(d)+d*dnorm(d))*sd^2
  
  # semisd
  semisd <- sd/sqrt(2)
  semimean <- -dnorm(0)*sd

  # VaR
  fq.alpha <- dnorm(q.alpha, mu, sd)

  # Nuisance Parameters for Performance Estimators
  
  # DSR
  dsr <- mu/semisd
  
  # ESratio
  ES.ratio <- mu/es.alpha
  
  # Omega
  upm1 <- lpm1+mu-c
  omega <- upm1/lpm1
  
  # RachR
  q.beta <- mu+qnorm(1-beta)*sd
  eg.beta <- mu+dnorm(qnorm(1-beta))/beta*sd
  rach.r <- eg.beta/es.alpha
  
  # robLoc
  psi.prime <- 0.9038

  # SoR_c
  sor.c <- mu/sqrt(lpm2)
  
  # SoR_mu
  sor.mu <- mu/semisd
  
  # SR 
  sr <- mu/sd   

  # VaRratio
  VaR.ratio <- -mu/q.alpha

  # Return list for nuisance parameters
  return(list(mu=mu, sd=sd, alpha=alpha, beta=beta, q.alpha=q.alpha, es.alpha=es.alpha, lpm1=lpm1, lpm2=lpm2, semisd=semisd, semimean=semimean, 
              fq.alpha=fq.alpha, dsr=dsr, ES.ratio=ES.ratio, upm1=upm1, omega=omega, c=c, q.beta=q.beta, eg.beta=eg.beta, rach.r=rach.r, psi.prime=psi.prime,
              sor.c=sor.c, sor.mu=sor.mu, sr=sr, VaR.ratio=VaR.ratio))
}