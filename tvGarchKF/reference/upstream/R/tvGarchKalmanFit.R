#' @title Fit the time-varying (Tv) parameters of the GARCH model (tv-Garch) by using the Kalman Filter method. The tv-parameters are determined by deterministic functions of either linear or non-linear type.
#' @description The tv-Garch(1,1) model, the parameters vary slowly over time according to linear or non-linear functions.
#' These parameters are denoted by \eqn{c(t)}, \eqn{\alpha(t)} and \eqn{\beta(t)} which correspond to the model
#' \eqn{\sigma_t = c(t) +\alpha(t) r^2_{t-1} +\beta(t)\sigma_{t-1}}.
#' @param series Time series.
#' @param c Vector containing coefficents of c.
#' @param alpha Vector containing coefficents of alpha.
#' @param beta Vector containing coefficents of beta.
#' @param type Vector of function type for c, alpha and beta.
#' @param exponentes Vector for exponenets in NoLineal.
#' @param trig Type of trigonometric function.
#' @param arg Value of argument for the trigonometric function.
#' @param predict Value for time to generate predict.
#' @param trace.log Variable to print names of coefficients.
#' @return Return fit values of omega, alpha and beta
#' @importFrom stats optim
#' @useDynLib tvGarchKF, .registration = TRUE
#' @export
#' @details
#' The types of functions for the tv-parameters are: linear, non-linear, trigonometric, and exponential.
#' For the case of the linear model, the tv-parameters follow the following structure:
#' \deqn{c(t) = c_0 + c_1u + c_2 u^2 + \ldots + c_p u^p,}
#' \deqn{\alpha(t) = a_0 + a_1u + a_2 u^2 + \ldots + a_p u^p,}
#' \deqn{\beta(t) = b_0 + b_1u + b_2 u^2 + \ldots + b_p u^p,}
#' where \eqn{u=t/T}, with \eqn{t=1, 2, \ldots, T}.
#' For the non-linear case, it is as follows:
#' \deqn{c(t) = c_0 + \sum_{j=1}^k c_j u_{c,j},}
#' \deqn{\alpha(t) = a_0 + \sum_{j=1}^k a_j u_{\alpha,j},}
#' \deqn{\beta(t) = b_0 + \sum_{j=1}^k b_j u_{\beta,j},}
#' where \eqn{k} it is positive value and \eqn{u_{c,j}}, \eqn{u_{\alpha,j}} and \eqn{u_{\beta,j}} are non linear function set.
#' For the trigonometric case, it is as follows:
#' \deqn{c(t) = c_0 + c_1 g(u),}
#' \deqn{\alpha(t) = a_0 + a_1 g(u),}
#' \deqn{\beta(t) = b_0 + b_1 g(u),}
#' where \eqn{g(u)} it is a trigonometric function, \eqn{cos} or \eqn{sin}.
#' @examples
#' ipsa<-diff(log(indipsa))
#' c <- c(0.05,0.05)
#' alpha <- c(0.05,0.05)
#' beta <- c(0.05,0.05)
#' type_fit <- c("trigonometric","trigonometric","trigonometric")
#' fit<-tvGarchKalmanFit(ipsa,c=c,alpha=alpha,beta=beta,type=type_fit,trig="cos",arg="3*(1-log(u))")

tvGarchKalmanFit <- function(series, c, alpha, beta, type=c("polynomial","NoLineal","trigonometric"), exponentes, trig, arg, predict = 0, trace.log=FALSE){
  checkInput(type)
  pars = c(c, alpha, beta)
  if(!missing(predict)){
    predict <- predict
  }else{
    predict <- 0
  }
  if(type[1]!="polynomial" | type[2]!="polynomial" | type[3]!="polynomial"){
    aux <- optim(par = pars, fn = tvGarchKalmanLoglike, series = series, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, method = c("BFGS"))
  }else{
    aux <- optim(par = pars, fn = tvGarchKalmanLoglike, series = series, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, method = c("CG"))
  }
  par <- round(aux$par,5)
  if(trace.log == T){
    for(i in 1:max(length(c),length(alpha),length(beta))){
      if(i<=length(c)){
        names(pars)[i] <- paste0("c",i)
      }
      if(i<=length(alpha)){
        names(pars)[i+length(c)] <- paste0("a",i)
      }
      if(i<=length(beta)){
        names(pars)[i+length(c)+length(alpha)] <- paste0("b",i)
      }
    }
  }
  suppressWarnings(return(par))
}

