#' Realized GARCH-Ito Model with Options
#' @description Estimate model parameters for the Realized GARCH-Ito Model with Options
#' @param RV Time series of daily realized volatilities.
#' @param JV Time series of daily jump variations,
#' @param NV Time series of daily volatilities estimated using option data
#' @param homogeneous Whether to assume homogeneous error in the linear regression model between conditional
#' volatility of the realized GARCH-Ito model and volatility estimated from the option data, default is TRUE.
#' @return Estimated parameter values and daily conditional volatilities:
#' \describe{
#' \item{coefficients}{parameter estimates of the realized GARCH-Ito model}
#' \item{sigma}{daily conditional volatility estimates of the realized GARCH-Ito model}
#' \item{pred}{one-step-ahead predicted volatility value}
#' }
#' @references Song, X., Kim, D., Yuan, H., Cui, X., Lu, Z., Zhou, Y., & Wang, Y. (2020).
#' Volatility Analysis with Realized GARCH-Ito Models. Journal of Econometrics, in press.
#'
#' @export

RealizedEst_Option=function(RV=RV, JV=NULL, NV=NULL, homogeneous=TRUE){

  if (length(NV) == 0 ){
    stop("estimated volatilities using option data are not available, use method RealizedEst() instead")
  }

  if (length(RV) != length(NV)){
    stop("realized volatility from high-frequency data and volatility from option data should be of equal length")
  }

  if (length(JV) == 0 ){

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

    if (homogeneous == TRUE ){

      # specify likelihood function
      likelihood.function.LHO=function(theta_tau){

        omega_g=theta_tau[1]; alpha_g=theta_tau[2]; gamma=theta_tau[3]
        a=theta_tau[4]; b=theta_tau[5]; sigma_e=theta_tau[6]

        hi=omega_g/(1-alpha_g-gamma)
        likelihood=log(hi)+RV[1]/hi+log(sigma_e^2)+(NV[1]-b-a*hi)^2/sigma_e^2
          for (i in 2:n){
            hi=hi.cond.function.NJ(theta_tau, hi, RV[i-1])
            likelihood=likelihood+log(hi)+RV[i]/hi+log(sigma_e^2)+(NV[i]-b-a*hi)^2/sigma_e^2
          }
        return(likelihood)
      }
      model=stats::lm(NV~RV)
      sigma_e=summary(model)$sigma; a=unname(model$coefficients[2]); b=unname(model$coefficients[1])
      # sink("file")
      thetaLHO=Rsolnp::solnp(c(mean(RV),0.4,0.3,a,b,sigma_e),fun=likelihood.function.LHO,ineqfun=ineqn.NJ,ineqLB=c(0),ineqUB=c(1),
                     LB=c(0,0,0,-10,-10,0),UB=c(1,1,1,10,10,10))$par
      # sink()
      names(thetaLHO)=c("omega_g", "alpha_g", "gamma", "a", "b", "sigma_e")

      hiLHO=rep(0,n+1); hiLHO[1]=thetaLHO[1]/(1-thetaLHO[2]-thetaLHO[3])
        for (i in 2:(n+1)){
          hiLHO[i]=hi.cond.function.NJ(thetaLHO, hiLHO[i-1], RV[i-1])
        }
      return(list("coefficients"=thetaLHO, "sigma"=hiLHO[1:n], "pred"=hiLHO[n+1]))

    } else {

      # specify likelihood function
      likelihood.function.LHOS=function(theta_tau){

        omega_g=theta_tau[1]; alpha_g=theta_tau[2]; gamma=theta_tau[3]
        a=theta_tau[4]; b=theta_tau[5]; sigma_e=theta_tau[6]; zeta=theta_tau[7]

        hi=omega_g/(1-alpha_g-gamma)
        likelihood=log(hi)+RV[1]/hi+log(sigma_e^2*hi^zeta)+(NV[1]-b-a*hi)^2/sigma_e^2/hi^zeta
          for (i in 2:n){
            hi=hi.cond.function.NJ(theta_tau, hi, RV[i-1])
            likelihood=likelihood+log(hi)+RV[i]/hi+log(sigma_e^2*hi^zeta)+(NV[i]-b-a*hi)^2/sigma_e^2/hi^zeta
          }
          return(likelihood)
        }
      model=stats::lm(NV~RV)
      sigma_e=summary(model)$sigma; a=unname(model$coefficients[2]); b=unname(model$coefficients[1])
      # sink("file")
      thetaLHOS=Rsolnp::solnp(c(mean(RV),0.4,0.3,a,b,sigma_e,sqrt(0.5)),fun=likelihood.function.LHOS,ineqfun=ineqn.NJ,ineqLB=c(0),ineqUB=c(1),
                      LB=c(0,0,0,-10,-10,0,0),UB=c(1,1,1,10,10,10,10))$par
      # sink()
      names(thetaLHOS)=c("omega_g", "alpha_g", "gamma", "a", "b", "sigma_e", "zeta")

      hiLHO=rep(0,n+1); hiLHO[1]=thetaLHOS[1]/(1-thetaLHOS[2]-thetaLHOS[3])
        for (i in 2:(n+1)){
          hiLHO[i]=hi.cond.function.NJ(thetaLHOS, hiLHO[i-1], RV[i-1])
        }
      return(list("coefficients"=thetaLHOS, "sigma"=hiLHO[1:n], "pred"=hiLHO[n+1]))
    }

  } else {

    if (length(RV) != length(JV)){
      stop("realized volatility and jump variation series should be of equal length")
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

    if ( homogeneous == TRUE){

      # specify likelihood function
      likelihood.function.LHO=function(theta_tau){

        omega_g=theta_tau[1]; alpha_g=theta_tau[2]; beta_g=theta_tau[3]; gamma=theta_tau[4]
        a=theta_tau[5]; b=theta_tau[6]; sigma_e=theta_tau[7]

        hi=(omega_g+beta_g*stats::median(JV))/(1-alpha_g-gamma)
        likelihood=log(hi)+RV[1]/hi+log(sigma_e^2)+(NV[1]-b-a*hi)^2/sigma_e^2

        for (i in 2:n){
          hi=hi.cond.function(theta_tau, hi, RV[i-1], JV[i-1])
          likelihood=likelihood+log(hi)+RV[i]/hi+log(sigma_e^2)+(NV[i]-b-a*hi)^2/sigma_e^2
        }
        return(likelihood)
      }

      model=stats::lm(NV~RV)
      sigma_e=summary(model)$sigma; a=unname(model$coefficients[2]); b=unname(model$coefficients[1])
      # sink("file")
      thetaLHO=Rsolnp::solnp(c(mean(RV),0.4,0.5,0.3,a,b,sigma_e),fun=likelihood.function.LHO,ineqfun=ineqn,ineqLB=c(0),ineqUB=c(1),
                     LB=c(0,0,0,0,-10,-10,0),UB=c(1,1,10,1,10,10,10))$par
      # sink()
      names(thetaLHO)=c("omega_g", "alpha_g", "beta_g", "gamma", "a", "b", "sigma_e")

      hiLHO=rep(0,n+1); hiLHO[1]=(thetaLHO[1]+thetaLHO[3]*stats::median(JV))/(1-thetaLHO[2]-thetaLHO[4])
        for (i in 2:(n+1)){
          hiLHO[i]=hi.cond.function(thetaLHO, hiLHO[i-1], RV[i-1], JV[i-1])
        }
      return(list("coefficients"=thetaLHO, "sigma"=hiLHO[1:n], "pred"=hiLHO[n+1]))

    } else {

      # specify likelihood function
      likelihood.function.LHOS=function(theta_tau){

        omega_g=theta_tau[1]; alpha_g=theta_tau[2]; beta_g=theta_tau[3]; gamma=theta_tau[4]
        a=theta_tau[5]; b=theta_tau[6]; sigma_e=theta_tau[7]; zeta=theta_tau[8]

        hi=(omega_g+beta_g*stats::median(JV))/(1-alpha_g-gamma)
        likelihood=log(hi)+RV[1]/hi+log(sigma_e^2*hi^zeta)+(NV[1]-b-a*hi)^2/sigma_e^2/hi^zeta
          for (i in 2:n){
            hi=hi.cond.function(theta_tau, hi, RV[i-1], JV[i-1])
            likelihood=likelihood+log(hi)+RV[i]/hi+log(sigma_e^2*hi^zeta)+(NV[i]-b-a*hi)^2/sigma_e^2/hi^zeta
          }
        return(likelihood)
      }

      model=stats::lm(NV~RV)
      sigma_e=summary(model)$sigma; a=unname(model$coefficients[2]); b=unname(model$coefficients[1])
      # sink("file")
      thetaLHOS=Rsolnp::solnp(c(mean(RV),0.4,0.5,0.3,a,b,sigma_e,sqrt(0.5)),fun=likelihood.function.LHOS,ineqfun=ineqn,ineqLB=c(0),ineqUB=c(1),
                      LB=c(0,0,0,0,-10,-10,0,0),UB=c(1,1,10,1,10,10,10,10))$par
      # sink()
      names(thetaLHOS)=c("omega_g", "alpha_g", "beta_g", "gamma", "a", "b", "sigma_e", "zeta")

      hiLHO=rep(0,n+1); hiLHO[1]=(thetaLHOS[1]+thetaLHOS[3]*stats::median(JV))/(1-thetaLHOS[2]-thetaLHOS[4])
        for (i in 2:(n+1)){
          hiLHO[i]=hi.cond.function(thetaLHOS, hiLHO[i-1], RV[i-1], JV[i-1])
        }
      return(list("coefficients"=thetaLHOS, "sigma"=hiLHO[1:n], "pred"=hiLHO[n+1]))

    }
  }
}
