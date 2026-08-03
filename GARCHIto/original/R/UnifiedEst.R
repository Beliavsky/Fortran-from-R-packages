#' Unified GARCH-Ito Models
#' @description Estimate model parameters for the Unified GARCH-Ito Model.
#' @param RV Time series of daily realized volatilities.
#' @param return Time series of daily log returns.
#' @return Estimated parameter values and daily conditional volatilities:
#' \describe{
#' \item{coefficients}{parameter estimates of the realized GARCH-Ito model}
#' \item{sigma}{daily conditional volatility estimates of the realized GARCH-Ito model}
#' \item{pred}{one-step-ahead predicted volatility value}
#' }
#' @examples
#' sample_data
#' UnifiedEst(sample_data$RV, sample_data$return)
#' @references Kim, D. & Wang, Y. (2016). Unified discrete-time and continuous-time models and statistical
#' inferences for merged low-frequency and high-frequency financial data. Journal of Econometrics. 194:220-230.
#' @export

UnifiedEst=function(RV=RV, return=return){

  n=length(RV)

  if(length(RV) != length(return)){
    stop("realized volatility and return series should be of equal length")
  }

  # conditional integrated volatility dynamics
  hi.cond.function.0=function(theta_tau, h, z){
    theta_tau[1]+theta_tau[3]*h+theta_tau[2]*z^2
  }

  # inequality for the model parameters
  ineqn.0=function(theta_tau){
    omega_g=theta_tau[1]; beta_g=theta_tau[2]; gamma=theta_tau[3]
    z1=beta_g+gamma
    return(c(z1))
  }

  # specify the propsoed likelihood function
  likelihood.function.LH.0=function(theta_tau){
    omega_g=theta_tau[1]; beta_g=theta_tau[2]; gamma=theta_tau[3]
    hi=theta_tau[1]/(1-theta_tau[2]-theta_tau[3])
    likelihood=log(hi)+RV[1]/hi
    for (i in 2:n){
      hi=hi.cond.function.0(theta_tau, hi, return[i-1])
      likelihood=likelihood+log(hi)+RV[i]/hi
    }
    return(likelihood)
  }

  # maximize likelihood function to obtain parameter estimates
  # sink("file")
  thetaLH0=Rsolnp::solnp(c(mean(RV),0.4,0.3),fun=likelihood.function.LH.0,ineqfun=ineqn.0,ineqLB=c(0),ineqUB=c(1),
                 LB=rep(0,3),UB=rep(1,3))$par
  # sink()
  names(thetaLH0)=c("omega_g", "beta_g", "gamma")

  hiLH0=rep(0,n+1); hiLH0[1]=thetaLH0[1]/(1-thetaLH0[2]-thetaLH0[3])
    for (i in 2:(n+1)){
      hiLH0[i]=hi.cond.function.0(thetaLH0, hiLH0[i-1], return[i-1])
    }
  return(list("coefficients"=thetaLH0, "sigma"=hiLH0[1:n], "pred"=hiLH0[n+1]))
}
