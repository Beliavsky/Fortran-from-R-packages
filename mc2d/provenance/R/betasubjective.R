#' @title The BetaSubjective Distribution
#' @name BetaSubjective
#' @keywords distribution
#' @description Density, distribution function, quantile function and random generation
#' for the "Beta Subjective" distribution
#' @details
#' The Subjective beta distribution specifies a [stats::dbeta()] distribution defined by the minimum, most likely (mode), mean
#' and maximum values and can be used for fitting data for a variable that is bounded to the interval \eqn{[min, max]}. 
#' The shape parameters are calculated from the mode value and mean parameters. It can also be used to represent 
#' uncertainty in subjective expert estimates.
#' 
#' Define
#' \deqn{mid=(min+max)/2}{mid=(min+max)/2}
#' \deqn{a_{1}=2*\frac{(mean-min)*(mid-mode)}{((mean-mode)*(max-min))}}{a_1=2*(mean-min)*(mid-mode)/((mean-mode)*(max-min))}
#' \deqn{a_{2}=a_{1}*\frac{(max-mean)}{(mean-min)}}{a_2=a_1*(max-mean)/(mean-min)}
#' 
#' The subject beta distribution is a [stats::dbeta()] distribution defined on the \eqn{[min, max]} domain 
#' with parameter \eqn{shape1 = a_{1}} and \eqn{shape2 = a_{2}}.  
#' 
#' # Hence, it has density
#' # \deqn{f(x)=(x-min)^{(a_{1}-1)}*(max-x)^{(a_{2}-1)} / (B(a_{1},a_{2})*(max-min)^{(a_{1}+a_{2}-1)})}{f(x)=(x-min)^(a_1-1)*(max-x)^(a_2-1)/(B(a_1,a_2)*(max-min)^(a_1+a_2-1))}
#' 
#' # The cumulative distribution function is
#' # \deqn{F(x)=B_{z}(a_{1},a_{2})/B(a_{1},a_{2})=I_{z}(a_{1},a_{2})}{F(x)=B_z(a_1,a_2)/B(a_1,a_2)=I_z(a_1,a_2)}
#' # where \eqn{z=(x-min)/(max-min)}. Here B is the beta function and \eqn{B_z} is the incomplete beta function.
#' 
#' The parameter restrictions are:
#' \deqn{min <= mode <= max}{min <= mode <= max}
#' \deqn{min <= mean <= max}{min <= mean <= max}
#' If \eqn{mode > mean} then \eqn{mode > mid}, else \eqn{mode < mid}.
#' 
#' @author Yu Chen 
#' 
#' @usage dbetasubj(x, 
#'   min,
#'   mode,
#'   mean,
#'   max, 
#'   log = FALSE)
#' @usage pbetasubj(q, 
#'   min,
#'   mode,
#'   mean,
#'   max, 
#'   lower.tail = TRUE,
#'   log.p = FALSE
#' )
#' @usage qbetasubj(p, 
#'   min,
#'   mode,
#'   mean,
#'   max, 
#'   lower.tail = TRUE, 
#'   log.p = FALSE
#' )
#' @usage rbetasubj(n, 
#'   min,
#'   mode,
#'   mean,
#'   max
#' )
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations. 
#' @param min continuous boundary parameter min < max
#' @param mode continuous parameter \eqn{min < mode < max} and \eqn{mode \ne mean}.
#' @param mean continuous parameter min < mean < max
#' @param max continuous boundary parameter
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' 
#' @examples 
#' curve(dbetasubj(x, min=0, mode=1, mean=2, max=5), from=-1,to=6) 
#' pbetasubj(q = seq(0,5,0.01), 0, 1, 2, 5)
#' qbetasubj(p = seq(0,1,0.01), 0, 1, 2, 5)
#' rbetasubj(n = 1e7, 0, 1, 2, 5)
#' @aliases dbetasubj
#' @aliases pbetasubj
#' @aliases qbetasubj
#' @aliases rbetasubj
#' @export
#' @encoding UTF-8
dbetasubj <- function(x,
                      min,
                      mode,
                      mean,
                      max,
                      log = FALSE)
{
  
  #parameter restriction
  #check min < mode < max
  if(any(min > mode)) stop("some min values are larger than mode values.")
  if(any(mode > max)) stop("some mode values are larger than max values.")
  
  #check min < mean < max
  if(any(min > mean)) stop("some min value are larger than mean values.")
  if(any(mean > max)) stop("some mean value are larger than max values.")  
  
  
  #check mode ?= mean
  if(any(mode == mean)) warning("mode is equal to mean.")
  
  mid <- (min + max)/2
  #check mode and mid
  if(any(mode>mean & mode<=mid)) stop("Impossible set of parameters (If mode is larger than mean, then mode should be larger than mid = (min + max)/2.")
  if(any(mode<mean & mode>=mid)) stop("Impossible set of parameters (If mode is smaller than mean, then mode should be smaller than mid = (min + max)/2.")
  
  a1 <- 2 * (mean-min) * (mid-mode) / (mean-mode) / (max-min)
  a2 <- a1 * (max-mean) / (mean-min)
  d <- (x-min)^(a1-1) * (max-x)^(a2-1) / beta(a1, a2) / (max-min)^(a1+a2-1)
  
  if(log) d <- log(d)
  
  return(d)
}

#' @name pbetasubj
#' @rdname BetaSubjective
#' @export
pbetasubj <- function(q,
                      min,
                      mode,
                      mean,
                      max,
                      lower.tail = TRUE,
                      log.p = FALSE)
{
  #parameter restriction
  #check min < mode < max
  if(any(min > mode)) stop("some min values are larger than mode values.")
  if(any(mode > max)) stop("some mode values are larger than max values.")
  
  #check min < mean < max
  if(any(min > mean)) stop("some min value are larger than mean values.")
  if(any(mean > max)) stop("some mean value are larger than max values.")  
  
  
  #check mode ?= mean
  if(any(mode == mean)) warning("mode is equal to mean.")
  
  mid <- (min + max)/2
  #check mode and mid
  if(any(mode>mean & mode<=mid)) stop("Impossible set of parameters (If mode is larger than mean, then mode should be larger than mid = (min + max)/2.")
  if(any(mode<mean & mode>=mid)) stop("Impossible set of parameters (If mode is smaller than mean, then mode should be smaller than mid = (min + max)/2.")
  
  a1 <- 2* (mean-min) * (mid-mode) / (mean-mode) / (max-min)
  a2 <- a1 * (max-mean) / (mean-min)
  
  z <- (q-min)/(max-min)
  p <- pbeta(z, a1, a2)
  
  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  
  return(p)
}

#' @name qbetasubj
#' @rdname BetaSubjective
#' @export 
qbetasubj <- function(p,
                 min,
                 mode,
                 mean,
                 max,
                 lower.tail = TRUE,
                 log.p = FALSE)
{
  #parameter restriction
  #check min < mode < max
  if(any(min > mode)) stop("some min values are larger than mode values.")
  if(any(mode > max)) stop("some mode values are larger than max values.")
  
  #check min < mean < max
  if(any(min > mean)) stop("some min value are larger than mean values.")
  if(any(mean > max)) stop("some mean value are larger than max values.")  
  
  
  #check mode ?= mean
  if(any(mode == mean)) warning("mode is equal to mean.")
  
  mid <- (min + max)/2
  #check mode and mid
  if(any(mode>mean & mode<=mid)) stop("Impossible set of parameters (If mode is larger than mean, then mode should be larger than mid = (min + max)/2.")
  if(any(mode<mean & mode>=mid)) stop("Impossible set of parameters (If mode is smaller than mean, then mode should be smaller than mid = (min + max)/2.")
  
  a1 <- 2* (mean-min) * (mid-mode) / (mean-mode) / (max-min)
  a2 <- a1 * (max-mean) / (mean-min)
  
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  z <- qbeta(p, a1, a2)
  q <- z * (max-min) + min 
  
  return(q)
}

#' @name rbetasubj
#' @rdname BetaSubjective
#' @export
rbetasubj <- function(n,
                 min,
                 mode,
                 mean,
                 max)
{

  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) cannot be negative in rbetasubj")
  
  #parameter restriction
  #check min < mode < max
  if(any(min > mode)) stop("some min values are larger than mode values.")
  if(any(mode > max)) stop("some mode values are larger than max values.")
  
  #check min < mean < max
  if(any(min > mean)) stop("some min value are larger than mean values.")
  if(any(mean > max)) stop("some mean value are larger than max values.")  
  
  
  #check mode ?= mean
  if(any(mode == mean)) warning("mode is equal to mean.")
  
  mid <- (min + max)/2
  #check mode and mid
  if(any(mode>mean & mode<=mid)) stop("Impossible set of parameters (If mode is larger than mean, then mode should be larger than mid = (min + max)/2.")
  if(any(mode<mean & mode>=mid)) stop("Impossible set of parameters (If mode is smaller than mean, then mode should be smaller than mid = (min + max)/2.")
  
  a1 <- 2* (mean-min) * (mid-mode) / (mean-mode) / (max-min)
  a2 <- a1 * (max-mean) / (mean-min)
  
  r <- (max-min) * rbeta(n, a1, a2) + min
  
  return(r)
}





