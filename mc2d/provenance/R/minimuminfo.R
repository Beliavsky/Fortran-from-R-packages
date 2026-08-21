#' @title Minimum Quantile Information Distribution
#' @name MinimumQuantileInformation
#' @keywords distribution
#' @description Density, distribution function, quantile function and random generation
#' for Minimum Quantile Information distribution.
#' @usage dmqi(x, 
#'   mqi, 
#'   mqi.quantile = c(0.05, 0.5, 0.95),
#'   realization = NULL, 
#'   k = 0.1, 
#'   intrinsic = NA,
#'   log = FALSE)
#' @usage pmqi(q, 
#'   mqi, 
#'   mqi.quantile = c(0.05, 0.5, 0.95),
#'   realization = NULL,
#'   k = 0.1,
#'   intrinsic = NA,
#'   lower.tail = TRUE,
#'   log.p = FALSE
#' )
#' @usage qmqi(p, 
#'   mqi, 
#'   mqi.quantile = c(0.05, 0.5, 0.95),
#'   realization = NULL, 
#'   k = 0.1, 
#'   intrinsic = NA,
#'   lower.tail = TRUE, 
#'   log.p = FALSE
#' )
#' @usage rmqi(n, 
#'   mqi, 
#'   mqi.quantile = c(0.05, 0.5, 0.95),
#'   realization = NULL, 
#'   k=0.1, 
#'   intrinsic = NA
#' )
#' @param x,q Vector of quantiles
#' @param p Vector of probabilities.
#' @param n Number of observations. 
#' @param mqi Minimum quantile information
#' @param mqi.quantile The quantile of `mqi`. It's a vector of length 3. Default is `c(0.05, 0.5, 0.95)`, 
#' that is the 5th, 50th and 95th.
#' @param realization Default is `NULL`. If not `NULL`, used to define `L` or `U` (see details).
#' @param k Overshot, default value is 0.1.
#' @param intrinsic Use to specify a prior bounds of the intrinsic range. Default = `NA`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are `P[X <= x]` otherwise, `P[X > x]`.
#' @param log,log.p Logical; if `TRUE`, probabilities `p` are given as `log(p)`.
#' @details
#' \eqn{p_1}, \eqn{p_2}, and \eqn{p_3} are percentiles of a distribution with \eqn{p_1 < p_2 < p_3}. 
#' The interval \eqn{[L,U]} is given with:
#' \deqn{L = x_{p_{1}}}{L = x_p_1}
#' \deqn{U = x_{p_{3}}}{U = x_p_3}
#' 
#' The support of minimum quantile information distribution is determined by the intrinsic range:
#' \deqn{[L^{*}, U^{*}] = [L - k \times (U - L), U + k \times (U - L)]}{[L^*, U^*] = [L - k*(U - L), U + k*(U - L)]}
#' where \eqn{k} denotes an overshoot and is chosen by the analyst (usually \eqn{k = 10\%}{k = 10\%}, which is the default value).
#' 
#' Given the three values of quantile, \eqn{x_{p_1}}, \eqn{x_{p_2}} and \eqn{x_{p_3}}, 
#' and define \eqn{p_0 = 0}, \eqn{p_4 = 1}, \eqn{x_{p_0} = L^{*}} and \eqn{x_{p_4} = U^{*}}
#' the minimum quantile information distribution is given by:
#' 
#' Probability density function 
#' \deqn{f(x) = \frac{p_{i}-p_{i-1}}{x_{p_{i}}-x_{p_{i-1}}} \text{ for } x_{p_{i-1}} \le x < x_{p_{i}},
#'  i = 1,\dots,4}{f(x)=(p_i-p_(i-1))/(x_p_i-x_p_(i-1)) for x_p_(i-1)\le x_p_i, i=1,\dots,4}
#' \deqn{f(x) = 0, \text{ otherwise}}{f(x) = 0, otherwise}
#' 
#' Cumulative distribution function
#' \deqn{F(x) = 0 \text{ for } x < x_{p_{0}}}{F(x) = 0 for x < x_p_0}
#' \deqn{F(x) = \frac{p_{i}-p_{i-1}}{x_{p_{i}}-x_{p_{i-1}}}*(x-x_{p_{i-1}})+p_{i-1} \text{ for } x_{p_{i-1}} \le x < x_{p_{i}}, i = 1,\dots,4}{F(x) = (p_i-p_(i-1))/(x_p_i-x_p_(i-1))*(x-x_p_(i-1))+p_(i-1) for x_p_(i-1)\le x_p_i, i=1,\dots,4}
#' 
#' \deqn{F(x) = 1 \text{ for } x_{p_{4}}\le x}{F(x) = 1 for x_p_(4) \le x}
#' 
#' 
#' 
#' This distribution is usually used for expert elicitation.
#' If experts have realization information, then the range \eqn{[L,U]} is given by:
#' \deqn{L = \min(x_{p_{1}}, realization)}{L = min(x_p_1, realization)}
#' \deqn{U = \max(x_{p_{3}}, realization)}{U = max(x_p_3, realization)}
#' 
#' For some questions, experts may have information for the intrinsic range and set a prior intrinsic range (\eqn{L^*} and \eqn{U^*}).
#' 
#' NOTE that the function is vectorized only for x, q, p, n. As a consequence, it can't be used 
#' for variable other parameters.
#' 
#' @author Yu Chen and Arie Havelaar
#' 
#' @references 
#' Hanea, A. M., & Nane, G. F. (2021). An in-depth perspective on the classical model. In International Series in Operations Research & Management Science (pp. 225–256). Springer International Publishing.
#' @examples 
#' 
#' curve(dmqi(x, mqi=c(40,50,60), intrinsic=c(0,100)), from=0, to=100, type = "l", xlab="x",ylab="pdf")
#' curve(pmqi(x, mqi=c(40,50,60), intrinsic=c(0,100)), from=0, to=100, type = "l", xlab="x",ylab="cdf")
#' rmqi(n = 10, mqi=c(555, 575, 586))
#' @aliases dmqi
#' @aliases pmqi
#' @aliases qmqi
#' @aliases rmqi
#' @export
dmqi <- function(x, 
                mqi, 
                mqi.quantile = c(0.05, 0.5, 0.95),
                realization = NULL, 
                k = 0.1, 
                intrinsic = NA, 
                log = FALSE)
{
  try(x <- as.matrix(x, nrow = 1))
  try(mqi <- as.matrix(mqi, nrow = 1))
  try(mqi.quantile <- as.vector(mqi.quantile))
  
  if(length(mqi)!=3) stop("Please input correct vector in mqi argument (length: 3)")
  if(length(mqi.quantile)!=3) stop("Please input correct vector in mqi.quantile argument (length: 3)")
  if(!(mqi[1] <= mqi[2] & mqi[2] <= mqi[3])) stop("Please input the mqi vector in increasing order")
  
  L <- min(min(mqi), realization)
  U <- max(max(mqi), realization)
  
  if(length(intrinsic)!=2 & sum(!is.na(intrinsic))!=0) stop("Please input correct intrinsic vector")
  
  if(is.na(intrinsic[1])) intrinsic[1] <- L - k * (U - L)
  if(is.na(intrinsic[2])) intrinsic[2] <- U + k * (U - L)
  
  # density in each interval
  
  
  f <- c(mqi.quantile[1], 
         (mqi.quantile[2]-mqi.quantile[1]), 
         (mqi.quantile[3]-mqi.quantile[2]), 
         1 - mqi.quantile[3]) / 
      c(mqi[1]-intrinsic[1], 
        mqi[2]-mqi[1],
        mqi[3]-mqi[2],
        intrinsic[2]-mqi[3])
  f <- c(0, f, 0)
  
  value.interval <- c(intrinsic[1], mqi, intrinsic[2]) 
  
  if(min(value.interval)!=intrinsic[1] | max(value.interval)!=intrinsic[2]) {
    stop("mqi is out of intrinsic range.")}

  index <- findInterval(x, value.interval)
  d <- f[index+1]
  d[x==intrinsic[2]] <- f[5]
  
  if(log) d <- log(d)
  return(d)
}

#' @name pmqi
#' @rdname MinimumQuantileInformation
#' @export
pmqi <- function(q, 
                mqi, 
                mqi.quantile = c(0.05, 0.5, 0.95),
                realization = NULL,
                k = 0.1,
                intrinsic = NA,
                lower.tail = TRUE,
                log.p = FALSE)
{
  try(q <- as.matrix(q, nrow = 1))
  try(mqi <- as.matrix(mqi, nrow = 1))
  try(mqi.quantile <- as.vector(mqi.quantile))
  
  if(length(mqi)!=3) stop("Please input correct vector in mqi argument")
  if(length(mqi.quantile)!=3) stop("Please input correct vector in mqi.quantile argument")
  if(!(mqi[1] <= mqi[2] & mqi[2] <= mqi[3])) stop("Please input the mqi vector in increasing order")
  
  L <- min(min(mqi), realization)
  U <- max(max(mqi), realization)
  
  if(length(intrinsic)!=2 & sum(!is.na(intrinsic))!=0) stop("Please input correct intrinsic vector")
  
  if(is.na(intrinsic[1])) intrinsic[1] <- L - k * (U - L)
  if(is.na(intrinsic[2])) intrinsic[2] <- U + k * (U - L)
  
  # density in each interval
  f <- c(0, 
         mqi.quantile[1], 
         (mqi.quantile[2]-mqi.quantile[1]), 
         (mqi.quantile[3]-mqi.quantile[2]),
         1 - mqi.quantile[3],
         0)
  
  value.interval <- c(intrinsic[1], mqi, intrinsic[2]) 
  
  prob.interval <- c(0, mqi.quantile, 1)
  
  index <- findInterval(q, value.interval)
  index[which(index==0)] <- 1
  index[which(index==5)] <- 4
  
  p <- f[index+1]/(value.interval[index+1]-value.interval[index]) * (q - value.interval[index]) + prob.interval[index]
  p <- as.numeric(p)
  
  p[which(p<0)] <- 0 #x < L* will output p<0
  p[which(p>1)] <- 1 #x > U* will output p>1
  
  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  
  return(p)
}

#' @name qmqi
#' @rdname MinimumQuantileInformation
#' @export 
qmqi <- function(p,
                mqi,
                mqi.quantile = c(0.05, 0.5, 0.95),
                realization=NULL,
                k = 0.1,
                intrinsic = NA,
                lower.tail = TRUE,
                log.p = FALSE
)
{
  try(p <- as.matrix(p, nrow = 1))
  try(mqi <- as.matrix(mqi, nrow = 1))
  try(mqi.quantile <- as.vector(mqi.quantile))
  
  if(length(mqi)!=3) stop("Please input correct vector in mqi argument")
  if(length(mqi.quantile)!=3) stop("Please input correct vector in mqi.quantile argument")
  if(!(mqi[1] <= mqi[2] & mqi[2] <= mqi[3])) stop("Please input the mqi vector in increasing order")
  
  L <- min(min(mqi), realization)
  U <- max(max(mqi), realization)
  
  if(length(intrinsic)!=2 & sum(!is.na(intrinsic))!=0) stop("Please input correct intrinsic vector")
  
  if(is.na(intrinsic[1])) intrinsic[1] <- L - k * (U - L)
  if(is.na(intrinsic[2])) intrinsic[2] <- U + k * (U - L)
  
  # density in each interval
  f <- c(0, 
         mqi.quantile[1], 
         (mqi.quantile[2]-mqi.quantile[1]), 
         (mqi.quantile[3]-mqi.quantile[2]),
         1 - mqi.quantile[3],
         0)
  
  value.interval <- c(intrinsic[1], mqi, intrinsic[2]) 
  
  prob.interval <- c(0, mqi.quantile, 1)
  
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  
  index <- findInterval(p, prob.interval)
  
  q <- rep(0, length(p))
  
  q <- (p - prob.interval[index]) * (value.interval[(index+1)] - value.interval[index]) /
    (prob.interval[(index+1)] - prob.interval[index]) + value.interval[index]
  
  q[p==1] <- intrinsic[2] #For p=1
  q <- as.numeric(q)
  
  return(q)
}

#' @name rmqi
#' @rdname MinimumQuantileInformation
#' @export
rmqi <- function(n,
                mqi,
                mqi.quantile = c(0.05, 0.5, 0.95),
                realization = NULL,
                k = 0.1,
                intrinsic = NA)
{
  
  try(mqi <- as.matrix(mqi, nrow = 1))
  try(mqi.quantile <- as.vector(mqi.quantile))

  
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) cannot be negative in rmqi")
  
  if(length(mqi)!=3) stop("Please input correct vector in mqi argument")
  if(length(mqi.quantile)!=3) stop("Please input correct vector in mqi.quantile argument")
  if(!(mqi[1] <= mqi[2] & mqi[2] <= mqi[3])) stop("Please input the mqi vector in increasing order")
  
  samples <- runif(n, 0, 1)
  
  r <- qmqi(samples, 
                 mqi = mqi, 
                 mqi.quantile = mqi.quantile,
                 realization = realization, 
                 k = k, 
                 intrinsic = intrinsic,
                 log.p = FALSE)
  
  return(r)
}

