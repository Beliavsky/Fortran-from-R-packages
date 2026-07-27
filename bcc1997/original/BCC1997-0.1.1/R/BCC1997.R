#' @title  Calculation of Option Prices Based on a Universal Solution
#' @description This is a function to calculate the prices of European options based on the universal solution provided by Bakshi, Cao and Chen (1997) <doi:10.1111/j.1540-6261.1997.tb02749.x>. This solution takes stochastic volatility, stochastic interest and random jumps into consideration. Please cite their work if this package is used.
#' @param kappav Speed of convergence on variance
#' @param kappar Speed of convergence on risk free rate
#' @param thetav Long-term variance
#' @param thetar Long-term risk free rate
#' @param sigmav Volatility of variance
#' @param sigmar Volatility of risk free rate
#' @param muj Jump size
#' @param sigmaj Volatility of jumps
#' @param rho Correlation between underlying price and variance
#' @param lambda Jump intensity
#' @param S0 Initial/Current underlying price
#' @param K Strike price
#' @param V0 Initial/Current variance
#' @param R0 Initial/Current risk free rate
#' @param t Time to maturity
#' @return Call: return the price of the European call oprion
#' @return Put: return the price of the European put oprion
#' @note Please notice each parameter has its "reasonable range". e.g. volatilities cannot be zero or smaller than zero, please input 0.0000001 when they are zero.
#' @examples BCC(kappav=0,kappar=0,thetav=0,thetar=0,sigmav=0.0000001,sigmar=0.0000001,muj=0,
#' @examples          sigmaj=0.0000001,rho=0,lambda=0,S0=100,K=100,V0=0.04,R0=0.01,t=1)
#' @examples BCC(kappav=0.5,kappar=0,thetav=0.025,thetar=0,sigmav=0.09,sigmar=0.0000001,muj=0,
#' @examples          sigmaj=0.0000001,rho=0.1,lambda=0,S0=100,K=100,V0=0.04,R0=0.01,t=1)
#' @export

BCC <- function(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t){
  f1 <- function(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda,S0,K,V0,R0,t,phi){
    er <- sqrt(kappar^2-2*sigmar^2*1i*phi)
    ev <- sqrt((kappav-(1+1i*phi)*rho*sigmav)^2-1i*phi*(1+1i*phi)*sigmav^2)
    Bt <- exp(-R0*t)
    a <- -thetar/(sigmar^2)*(2*log(1-(er-kappar)*(1-exp(-er*t))/(2*er))+(er-kappar)*t)
    b <- -thetav/(sigmav^2)*2*log(1-(ev-kappav+(1+1i*phi)*rho*sigmav)*(1-exp(-ev*t))/(2*ev))
    c <- -thetav/(sigmav^2)*(ev-kappav+(1+1i*phi)*rho*sigmav)*t
    d <- 1i*phi*log(S0)
    e <- 2*1i*phi*(1-exp(-er*t))/(2*er-(er-kappar)*(1-exp(-er*t)))*R0
    f <- lambda*(1+muj)*t*((1+muj)^(1i*phi)*exp(1i*phi/2*(1+1i*phi)*sigmaj^2)-1)
    g <- -lambda*1i*phi*muj*t
    h <- V0*1i*phi*(1+1i*phi)*(1-exp(-ev*t))/(2*ev-(ev-kappav+(1+1i*phi)*rho*sigmav)*(1-exp(-ev*t)))
    return(exp(a+b+c+d+e+f+g+h))
  }

  f2 <- function(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda,S0,K,V0,R0,t,phi){
    er <- sqrt(kappar^2-2*sigmar^2*(1i*phi-1))
    ev <- sqrt((kappav-1i*phi*rho*sigmav)^2-1i*phi*(1i*phi-1)*sigmav^2)
    Bt <- exp(-R0*t)
    a <- -thetar/(sigmar^2)*(2*log(1-(er-kappar)*(1-exp(-er*t))/(2*er))+(er-kappar)*t)
    d <- 1i*phi*(log(S0))-log(Bt)
    b <- -thetav/(sigmav^2)*2*log(1-(ev-kappav+1i*phi*rho*sigmav)*(1-exp(-ev*t))/(2*ev))
    c <- -thetav/(sigmav^2)*(ev-kappav+1i*phi*rho*sigmav)*t
    g <- -lambda*1i*phi*muj*t
    e <- 2*(1i*phi-1)*(1-exp(-er*t))/(2*er-(er-kappar)*(1-exp(-er*t)))*R0
    f <- lambda*t*((1+muj)^(1i*phi)*exp(1i*phi/2*(1i*phi-1)*sigmaj^2)-1)
    h <- V0*1i*phi*(1i*phi-1)*(1-exp(-ev*t))/(2*ev-(ev-kappav+1i*phi*rho*sigmav)*(1-exp(-ev*t)))
    return(exp(a+b+c+d+e+f+g+h))
  }

  Pi1 <- function(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t){
    integrand <- function(phi){Re(exp(-1i*phi*log(K))*f1(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda,S0,K,V0,R0,t,phi)/(1i*phi))}
    vvpi <- 1/2+1/pi*integrate(integrand,lower=0,upper=Inf,subdivisions = 10000L)$value
    return(vvpi)
  }

  Pi2 <- function(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t){
    integrand <- function(phi){Re(exp(-1i*phi*log(K))*f2(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda,S0,K,V0,R0,t,phi)/(1i*phi))}
    vvpi <- 1/2+1/pi*integrate(integrand,lower=0,upper=Inf,subdivisions = 10000L)$value
    return(vvpi)
  }

  Bt <- exp(-R0*t)
  call <- S0*Pi1(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t)-K*Bt*Pi2(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t)
  put <- -S0*(1-Pi1(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t))+K*Bt*(1-Pi2(kappav,kappar,thetav,thetar,sigmav,sigmar,muj,sigmaj,rho,lambda, S0,K,V0,R0,t))
  result <- list(call=call,put=put)
  return(result)
}

