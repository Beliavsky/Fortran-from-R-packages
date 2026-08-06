#' @title Models tv-Garch Filter Kalman print outputs.
#' @description This function is designed to print the outputs of the tv-Garch model, which include the returns,
#'    conditional variance, log-likelihood value, and mean squared error (MSE).
#' @param x Vector of coefficents to fit.
#' @param series Time series.
#' @param c Vector containing coefficents of c.
#' @param alpha Vector containing coefficents of alpha.
#' @param beta Vector containing coefficents of beta.
#' @param nsample Value of time series length.
#' @param type Vector of function type for c, alpha and beta.
#' @param exponentes Vector for exponenets in NoLineal.
#' @param trig Type of trigonometric function.
#' @param arg Value of argument for the trigonometric function.
#' @param predict Value for time to generate predict.
#' @param trace.log Variable to print names of coefficients.
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item \code{X}: State vector of Kalman equations.
#'   \item \code{Fm}: Value of MSE
#'   \item \code{sigma}: Conditional variance.
#'   \item \code{loglike}: Value of the loglike.
#' }
#' @importFrom stats na.exclude
#' @export
#' @examples
#' data(ipsa)
#' ipsa<-diff(log(indipsa))
#' c<-c(0.05,0.05)
#' alpha<-c(0.05,0.05)
#' beta<-c(0.05,0.05)
#' type_fit<-c("trigonometric","trigonometric","trigonometric")
#' fit<-tvGarchKalmanFit(ipsa,c=c,alpha=alpha,beta=beta,type=type_fit,trig="cos",arg="3*(1-log(u))")
#' arg_model<-"3*(1-log(u))"
#' model<-tvGarchKalmanPrint(fit,ipsa,c=c,alpha=alpha,beta=beta,type=type_fit,trig="cos",arg=arg_model)
#' plot(ipsa,ylab="",xlim=c(2000,2015))
#' lines(ts(model$sigma, star=2000, freq=225), col="red", lwd=2)
#' lines(ts(model$sigma*(-1), star=2000, freq=225), col="red", lwd=2)
tvGarchKalmanPrint <- function(x, series, c, alpha, beta, nsample=length(series), type=c("polynomial","NoLineal","trigonometric"), exponentes, trig, arg, trace.log=FALSE, predict){

  pred<-rep(NA,predict)
  series<-c(series,pred)

  nsample=length(series)
  u <- (1:nsample)/nsample

  aux1 <- length(c)
  aux2 <- length(alpha)
  aux3 <- length(beta)
  if (aux1<=1){
    c = x[1]
  }else{
    c = c()
    for(i in 1:aux1){
      c[i] = x[i]
    }
  }
  if (aux2<=1){
    alpha = x[1+aux1]
  }else{
    alpha = c()
    for(i in 1:aux2){
      alpha[i] = x[i+aux1]
    }
  }
  if (aux3<=1){
    beta = x[1+aux1+aux2]
  }else{
    beta = c()
    for(i in 1:aux3){
      beta[i] = x[i+aux2+aux1]
    }
  }

  auxA <- list(c=c,alpha=alpha,beta=beta)
  auxB <- list()
  for(i in 1:length(type)){
    if(type[i] == "polynomial"){
      auxC <- polynomial(nsample, coeficientes = unlist(auxA[i]))
    }
    if(type[i] == "NoLineal"){
      auxC <- NoLineal(nsample, coeficientes = unlist(auxA[i]), exponentes = exponentes)
    }
    if(type[i] == "trigonometric"){
      auxC <- trigonometric(nsample, coeficientes = unlist(auxA[i]), trig = trig, arg = arg)
    }
    auxB <- c(auxB,list(auxC))
  }
  c <- unlist(auxB[1])
  a <- unlist(auxB[2])
  b <- unlist(auxB[3])


  for(i in 1:nsample){
    if(any(is.nan(a[i])) || any(is.nan(b[i])) || any(is.nan(c))) {
      stop("Se produjo un error al estimar los par?metros")
    }
  }


  for(i in 1:nsample){
    if(trace.log) {
      cat("a[i]:",a[i],"|b[i]:",b[i],"|c:",c,"\n")
    }
  }

  for(i in 1:nsample){
    if ((c[i]<=0) & (a[i]<=0) & (b[i]<=0) & ((b[i]+a[i])>=1) ) {
      return(100000000.)
    }
  }

  y<-NULL
  for(i in 1:nsample){
    y[i] <- ifelse (is.na(series[i]), 100000000., series[i]^2-c[i]/(1-a[i]-b[i]))
  }


  P <-rep(0, nsample) # Varianza del error de estimaci?n
  X <-rep(0, nsample)   # Representa el estado (alfa)
  F.m <-rep(0, nsample)
  K<-rep(0, nsample)

  # Inicializo P_1
  P[1]<-a[1]^2/(1-(a[1] + b[1])^2)

  if (aux1 > 1){
    resultado<-.C("tvGarch1ab2", as.integer(nsample), as.double(a), as.double(b), as.double(c), as.double(P), as.double(X), as.double(F.m), as.double(K), as.double(y))
  }else{
    resultado<-.C("tvGarch1ab", as.integer(nsample), as.double(a), as.double(b), as.double(c), as.double(P), as.double(X), as.double(F.m), as.double(K), as.double(y))
  }

  X=resultado[[6]]
  F.m=resultado[[7]]
  K=resultado[[8]]

  ####

  e.sig2 = X + (c/(1-a-b))

  if( min(e.sig2) <= 0 | min(F.m) <= 0 ){
    loglik=100000000
    }
  else{
    loglik <- sum(na.exclude((sqrt(F.m)*series^2) / e.sig2 + log(e.sig2) - 0.5*log(F.m) ))
    loglik
    }


  OUT = NULL
  if (all(e.sig2>=0)){
    sigma<-sqrt(e.sig2)
    OUT$sigma = sigma
  }else{print("Hay valores negativos en sigma.")}
  OUT$X=X
  OUT$Fm=F.m
  OUT$loglik = loglik
  OUT

}
