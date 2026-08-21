#'Autocorrelation function of the exponential trawl function
#'
#'This function computes the autocorrelation function associated with the
#'exponential trawl function.
#'@name acf_Exp
#'@param x The argument (lag) at which the autocorrelation function associated
#'  with the exponential trawl function will be evaluated
#'@param lambda parameter in the exponential trawl
#'@return The autocorrelation function of the exponential trawl function
#'  evaluated at x
#'@details The trawl function is parametrised by the parameter \eqn{\lambda > 0}
#'  as follows:
#'\deqn{g(x) = e^{\lambda x},  \mbox{ for } x \le 0.}
#'Its autocorrelation function is given by:
#'\deqn{r(x) = e^{-\lambda x},  \mbox{ for } x \ge 0.}
#'@examples
#'acf_Exp(1,0.1)
#'@export
acf_Exp <- function(x,lambda){exp(-lambda*x)}


#'Autocorrelation function of the double exponential trawl function
#'
#'This function computes the autocorrelation function associated with the double
#'exponential trawl function.
#'@name acf_DExp
#'@param x The argument (lag) at which the autocorrelation function associated
#'  with the double exponential trawl function will be evaluated
#'@param w parameter in the double exponential trawl
#'@param lambda1 parameter in the double exponential trawl
#'@param lambda2 parameter in the double exponential trawl
#'@return The autocorrelation function of the double exponential trawl function
#'  evaluated at x
#'@details The trawl function is parametrised by parameters \eqn{0\le w\le 1}
#'  and \eqn{\lambda_1, \lambda_2 > 0}  as follows: \deqn{g(x) = w e^{\lambda_1
#'  x}+(1-w) e^{\lambda_2 x},  \mbox{ for }  x \le 0.} Its autocorrelation
#'  function is given by: \deqn{r(x) = (w e^{-\lambda_1 x}/\lambda_1+(1-w)
#'  e^{-\lambda_2 x}/\lambda_2)/c,  \mbox{ for } x \ge 0,} where \deqn{c =
#'  w/\lambda_1+(1-w)/\lambda_2.}
#'@examples
#'acf_DExp(1,0.3,0.1,2)
#'@export
acf_DExp <- function(x,w,lambda1,lambda2){(w*exp(-lambda1*x)/lambda1+(1-w)*exp(-lambda2*x)/lambda2)/(w/lambda1+(1-w)/lambda2)}

#'Autocorrelation function of the supIG trawl function
#'
#'This function computes the autocorrelation function associated with the supIG
#'trawl function.
#'@name acf_supIG
#'@param x The argument (lag) at which the autocorrelation function associated
#'  with the supIG trawl function will be evaluated
#'@param delta parameter in the supIG trawl
#'@param gamma parameter in the supIG trawl
#'@return The autocorrelation function of the supIG trawl function evaluated at
#'  x
#'@details The trawl function is parametrised by the two parameters \eqn{\delta
#'  \geq 0} and \eqn{\gamma \geq 0} as follows: \deqn{g(x) =
#'  (1-2x\gamma^{-2})^{-1/2}\exp(\delta \gamma(1-(1-2x\gamma^{-2})^{1/2})),
#'  \mbox{ for } x \le 0.} It is assumed that \eqn{\delta} and \eqn{\gamma} are
#'  not simultaneously equal to zero. Its autocorrelation function is given by:
#'  \deqn{r(x) = \exp(\delta\gamma (1-\sqrt{1+2 x/\gamma^2})),  \mbox{ for } x
#'  \ge 0.}
#'@examples
#'acf_supIG(1,0.3,0.1)
#'@export
acf_supIG <- function(x,delta,gamma){exp(delta*gamma*(1-sqrt(1+2*x/gamma^2)))}

#'Autocorrelation function of the long memory trawl function
#'
#'This function computes the autocorrelation function associated with the long
#'memory trawl function.
#'@name acf_LM
#'@param x The argument (lag) at which the autocorrelation function associated
#'  with the long memory trawl function will be evaluated
#'@param alpha parameter in the long memory trawl
#'@param H parameter in the long memory trawl
#'@return The autocorrelation function of the long memory trawl function
#'  evaluated at x
#'@details The trawl function is parametrised by the two parameters \eqn{H> 1}
#'  and \eqn{\alpha > 0} as follows: \deqn{g(x) = (1-x/\alpha)^{-H},  \mbox{ for
#'  }  x \le 0.} Its autocorrelation function is given by
#'  \deqn{r(x)=(1+x/\alpha)^{(1-H)}, \mbox{ for } x \ge 0.}
#'@examples
#'acf_LM(1,0.3,1.5)
#'@export
acf_LM <- function(x,alpha, H){(1+x/alpha)^(1-H)}
