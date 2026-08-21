#'Fits an exponential trawl function to equidistant time series data
#'@name fit_Exptrawl
#'@param x vector of equidistant time series data
#'@param Delta interval length of the time grid used in the time series, the
#'  default is 1
#'@param plotacf binary variable specifying whether or not the empirical and
#'  fitted autocorrelation function should be plotted
#'@param lags number of lags to be used in the plot of the autocorrelation
#'  function
#'@return lambda: the memory parameter \eqn{\lambda} in the exponential trawl
#'@return LM: the Lebesgue measure of the trawl set associated with the
#'  exponential trawl,  i.e. \eqn{1/\lambda}.
#'@details The trawl function is parametrised by the parameter \eqn{\lambda > 0}
#'  as follows: \deqn{g(x) = e^{\lambda x},  \mbox{ for }  x \le 0.} The
#'  Lebesgue measure of the corresponding trawl set is given by \eqn{1/\lambda}.
#'@export
fit_Exptrawl <- function(x,Delta=1, plotacf=FALSE,lags=100){
  my_acf <- TSA::acf(x,plot=F)
  lambda <- -base::log(my_acf$acf[1])/Delta

  if(plotacf){

    tt <- (1:lags)
    TSA::acf(x,lag.max=lags,main = "", ylab="ACF", xlab="Lags",drop.lag.0 = TRUE)
    graphics::lines(tt, acf_Exp((tt*Delta), lambda), lty =1,col=2, lwd=2)
  }

  return (list("lambda"=lambda,"LM"=1/lambda))
}






#'Fits a supIG trawl function to equidistant univariate time series data
#'@name fit_supIGtrawl
#'@param x vector of equidistant time series data
#'@param Delta interval length of the time grid used in the time series, the
#'  default is 1
#'@param GMMlag lag length used in the GMM estimation, the default is 5
#'@param plotacf binary variable specifying whether or not the empirical and
#'  fitted autocorrelation function should be plotted
#'@param lags number of lags to be used in the plot of the autocorrelation
#'  function
#'@return delta: parameter in the supIG trawl
#'@return gamma: parameter in the supIG trawl
#'@return LM: The Lebesgue measure of the trawl set associated with the supIG
#'  trawl
#'@details The trawl function is parametrised by the two parameters \eqn{\delta
#'  \geq 0} and \eqn{\gamma \geq 0} as follows: \deqn{g(x) =
#'  (1-2x\gamma^{-2})^{-1/2}\exp(\delta \gamma(1-(1-2x\gamma^{-2})^{1/2})),
#'  \mbox{ for }  x \le 0.} It is assumed that \eqn{\delta} and \eqn{\gamma} are
#'  not simultaneously equal to zero. The Lebesgue measure of the corresponding
#'  trawl set is given by \eqn{\gamma/\delta}.
#'@export
fit_supIGtrawl <- function(x,Delta=1,GMMlag=5,plotacf=FALSE,lags=100){
  my_acf <- TSA::acf(x,plot=F)

  fit_supIGtrawl_foroptim <- function(y){

    delta <- y[1]
    gamma <-y[2]

    lag <- GMMlag
    lss <- 0

    for(i in 1:lag)
    {
      Cor <- acf_supIG(Delta*i,delta,gamma)
      lss <- lss+(my_acf$acf[i] - Cor)^2
    }
    lss
  }

  o <- DEoptim::DEoptim(fit_supIGtrawl_foroptim,c(0,0),c(100,100),control=DEoptim::DEoptim.control(itermax = 1000, trace = FALSE))

  if(plotacf){

    tt <- (1:lags)
    TSA::acf(x,lag.max=lags,main = "", ylab="ACF", xlab="Lags",drop.lag.0 = TRUE)
    graphics::lines(tt, acf_supIG(tt*Delta,o$optim$bestmem[1],o$optim$bestmem[2]), lty =1,col=2, lwd=2)
  }

  return(list("delta"=as.numeric(o$optim$bestmem[1]), "gamma"=as.numeric(o$optim$bestmem[2]), "LM"=as.numeric(o$optim$bestmem[2])/as.numeric(o$optim$bestmem[1])))
}





#'Fits a long memory trawl function to equidistant univariate time series data
#'@name fit_LMtrawl
#'@param x vector of equidistant time series data
#'@param Delta interval length of the time grid used in the time series, the
#'  default is 1
#'@param GMMlag lag length used in the GMM estimation, the default is 5
#'@param plotacf binary variable specifying whether or not the empirical and
#'  fitted autocorrelation function should be plotted
#'@param lags number of lags to be used in the plot of the autocorrelation
#'  function
#'@return alpha: parameter in the long memory trawl
#'@return H: parameter in the long memory trawl
#'@return LM: The Lebesgue measure of the trawl set associated with the
#' long memory  trawl
#'@details The trawl function is parametrised by the two parameters \eqn{H> 1}
#'  and \eqn{\alpha > 0} as follows: \deqn{g(x) = (1-x/\alpha)^{-H},\mbox{ for }
#'  x \le  0.} The Lebesgue measure of the corresponding trawl set is given by
#'  \eqn{\alpha/(1-H)}.
#'@export
fit_LMtrawl <- function(x,Delta=1, GMMlag=5, plotacf=FALSE,lags=100){
  my_acf <- TSA::acf(x,plot=F)

  fit_LMtrawl_foroptim <- function(y){

    alpha <- y[1]
    Hurst <-y[2]

    lag <- GMMlag
    lss <- 0

    for(i in 1:lag)
    {
      Cor <- acf_LM(Delta*i,alpha, Hurst)
      lss <- lss+(my_acf$acf[i] - Cor)^2
    }

    lss

  }

  o <- DEoptim::DEoptim(fit_LMtrawl_foroptim,c(0,0),c(100,100),control=DEoptim::DEoptim.control(itermax = 1000, trace = FALSE))

  if(plotacf){

    tt <- (1:lags)
    TSA::acf(x,lag.max=lags,main = "", ylab="ACF", xlab="Lags",drop.lag.0 = TRUE)
    graphics::lines(tt, acf_LM(tt*Delta,o$optim$bestmem[1],o$optim$bestmem[2]), lty =1,col=2, lwd=2)
  }

  return(list("alpha"=as.numeric(o$optim$bestmem[1]), "H"=as.numeric(o$optim$bestmem[2]), "LM"=as.numeric(o$optim$bestmem[1])/(as.numeric(o$optim$bestmem[2])-1)))
}




#'Fits the trawl function consisting of the weighted sum of two exponential functions
#'@name fit_DExptrawl
#'@param x vector of equidistant time series data
#'@param Delta interval length of the time grid used in the time series, the
#'  default is 1
#'@param GMMlag lag length used in the GMM estimation, the default is 5
#'@param plotacf binary variable specifying whether or not the empirical and
#'  fitted autocorrelation function should be plotted
#'@param lags number of lags to be used in the plot of the autocorrelation
#'  function
#'@return w: the weight parameter (restricted to be in [0,0.5] for
#'  identifiability reasons)
#'@return lambda1: the first memory parameter (denoted by \eqn{\lambda_1} above)
#'@return lambda2: the second memory parameter (denoted by \eqn{\lambda_2}
#'  above)
#'@return LM: The Lebesgue measure of the trawl set associated with the double
#'  exponential trawl
#'@details The trawl function is parametrised by the three parameters \eqn{0\leq
#'  w \leq 1} and \eqn{\lambda_1,\lambda_2 > 0} as follows: \deqn{g(x) =
#'  we^{\lambda_1 x}+(1-w)e^{\lambda_2 x},  \mbox{ for }  x \le 0.} The Lebesgue measure
#'  of the corresponding trawl set is given by
#'  \eqn{w/\lambda_1+(1-w)/\lambda_2}.
#'@export
fit_DExptrawl <- function(x,Delta=1, GMMlag=5, plotacf=FALSE,lags=100){
  my_acf <- TSA::acf(x,plot=F)

  fit_DExptrawl_foroptim <- function(y){
    w <- y[1]
    memory1 <-y[2]
    memory2 <-y[3]

    Delta <- 1
    my_acf <- TSA::acf(x,plot=F)

    lag <- GMMlag
    lss <- 0

    for(i in 1:lag)
    {
      Cor <- acf_DExp(Delta*i,w,memory1,memory2)
      lss <- lss+(my_acf$acf[i] - Cor)^2
    }

    lss
  }

  o <- DEoptim::DEoptim(fit_DExptrawl_foroptim,c(0,0,0),c(0.5,100,100),control=DEoptim::DEoptim.control(itermax = 1000, trace = FALSE))

  if(plotacf){

    tt <- (1:lags)
    TSA::acf(x,lag.max=lags,main = "", ylab="ACF", xlab="Lags",drop.lag.0 = TRUE)
    graphics::lines(tt, acf_DExp(tt*Delta,o$optim$bestmem[1],o$optim$bestmem[2],o$optim$bestmem[3]), lty =1,col=2, lwd=2)
  }

  return(list("w"=as.numeric(o$optim$bestmem[1]), "lambda_1"=as.numeric(o$optim$bestmem[2]), "lambda_2"=as.numeric(o$optim$bestmem[3]),"LM"=as.numeric(o$optim$bestmem[1])/as.numeric(o$optim$bestmem[2])+(1-as.numeric(o$optim$bestmem[1]))/as.numeric(o$optim$bestmem[3])))
}

