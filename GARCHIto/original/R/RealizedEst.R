#' Realized GARCH-Ito Model
#' @description Estimate model parameters for the Realized GARCH-Ito Model
#' @param RV Time series of daily realized volatilities.
#' @param JV Time series of daily jump variations,
#' @return Estimated parameter values and daily conditional volatilities:
#' \describe{
#' \item{coefficients}{parameter estimates of the realized GARCH-Ito model}
#' \item{sigma}{daily conditional volatility estimates of the realized GARCH-Ito model}
#' \item{pred}{one-step-ahead predicted volatility value}
#' }
#' @examples
#' sample_data
#' RealizedEst(sample_data$RV)
#' RealizedEst(sample_data$BPV, sample_data$JV)
#' @references Song, X., Kim, D., Yuan, H., Cui, X., Lu, Z., Zhou, Y., & Wang, Y. (2020).
#' Volatility Analysis with Realized GARCH-Ito Models. Journal of Econometrics, in press.
#' @export

RealizedEst=function(RV=RV, JV=NULL){

  if (length(JV) == 0){

    n=length(RV)

    # conditional integrated volatility dynamics
    hi.cond.function.NJ=function(theta_tau, h, iv){
      theta_tau[1]+theta_tau[3]*h+theta_tau[2]*iv
    }

    # constraint for GARCH model parameters
    ineqn.NJ=function(theta_tau){
      omega_g=theta_tau[1]; alpha_g=theta_tau[2]; gamma=theta_tau[3]
      z1=alpha_g+gamma
      return(c(z1))
    }

    # specify likelihood function
    likelihood.function.LH.NJ=function(theta_tau){
      omega_g=theta_tau[1]; alpha_g=theta_tau[2]; gamma=theta_tau[3]
      hi=omega_g/(1-alpha_g-gamma)
      likelihood=log(hi)+RV[1]/hi
        for (i in 2:n){
          hi=hi.cond.function.NJ(theta_tau, hi, RV[i-1])
          likelihood=likelihood+log(hi)+RV[i]/hi
        }
      return(likelihood)
    }

    # maximize likelihood function to obtain parameter estimates
    # sink("file")
    thetaLH=Rsolnp::solnp(c(mean(RV),0.45,0.3),fun=likelihood.function.LH.NJ,ineqfun=ineqn.NJ,ineqLB=c(0),ineqUB=c(1),
                         LB=rep(0,3),UB=rep(1,3))$par
    # sink()
    names(thetaLH)=c("omega_g", "alpha_g", "gamma")

    hiLH=rep(0,n+1); hiLH[1]=thetaLH[1]/(1-thetaLH[2]-thetaLH[3])
      for (i in 2:(n+1)){
        hiLH[i]=hi.cond.function.NJ(thetaLH, hiLH[i-1], RV[i-1])
      }
    return(list("coefficients"=thetaLH, "sigma"=hiLH[1:n], "pred"=hiLH[n+1]))

  } else {

    if (length(RV) != length(JV)){
      stop("realized volatility and jump variation series should be of equal length")
    }

    if (length(which(RV<0)) != 0){
      stop("realized volatility should be non-negative")
    }

    if (length(which(JV<0)) != 0){
      stop("jump variation should be non-negative")
    }

    n=length(RV)

    # conditional integrated volatility dynamics
    hi.cond.function=function(theta_tau, h, iv, jv){
      theta_tau[1]+theta_tau[4]*h+theta_tau[2]*iv+theta_tau[3]*jv
    }

    # constraint for GARCH model parameters
    ineqn=function(theta_tau){
      omega_g=theta_tau[1]; alpha_g=theta_tau[2]; beta_g=theta_tau[3]; gamma=theta_tau[4]
      z1=alpha_g+gamma
      return(c(z1))
    }

    # specify likelihood function
    likelihood.function.LH=function(theta_tau){

      omega_g=theta_tau[1]; alpha_g=theta_tau[2]; beta_g=theta_tau[3]; gamma=theta_tau[4]
      hi=(omega_g+beta_g*stats::median(JV))/(1-alpha_g-gamma)
      likelihood=log(hi)+RV[1]/hi
        for (i in 2:n){
          hi=hi.cond.function(theta_tau, hi, RV[i-1], JV[i-1])
          likelihood=likelihood+log(hi)+RV[i]/hi
        }
        return(likelihood)
      }

    # maximize likelihood function to obtain parameter estimates
    # sink("file")
    thetaLH=Rsolnp::solnp(c(mean(RV),0.45,0.5,0.3),fun=likelihood.function.LH,ineqfun=ineqn,ineqLB=c(0),ineqUB=c(1),
                          LB=c(0,0,0,0),UB=c(1,1,10,1))$par
    # sink()
    names(thetaLH)=c("omega_g", "alpha_g", "beta_g", "gamma")

    hiLH=rep(0,n+1); hiLH[1]=(thetaLH[1]+thetaLH[3]*stats::median(JV))/(1-thetaLH[2]-thetaLH[4])
      for (i in 2:(n+1)){
        hiLH[i]=hi.cond.function(thetaLH, hiLH[i-1], RV[i-1], JV[i-1])
      }
    return(list("coefficients"=thetaLH, "sigma"=hiLH[1:n], "pred"=hiLH[n+1]))
  }

}
