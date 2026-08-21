#' The Generalised Beta Distribution
#'
#' Density, distribution function, quantile function and random generation for the Beta
#' distribution defined on the \samp{[min, max]} domain with parameters \samp{shape1}
#' and \samp{shape2} (and optional non-centrality parameter \samp{ncp}).
#'
#' \deqn{x \sim betagen(shape1, shape2, min, max, ncp)}{x ~ betagen(shape1, shape2, min, max, ncp)}
#' if
#' \deqn{\frac{x-min}{max-min}\sim beta(shape1,shape2,ncp)}{(x-min)/(max-min)~beta(shape1,shape2,ncp)}
#' These functions use the \code{\link{Beta}} distribution functions after correct parameterization.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \samp{length(n) > 1}, the length is taken to be the number required.
#' @param shape1,shape2 positive parameters of the Beta distribution.
#' @param min vector of minima.
#' @param max vector of maxima.
#' @param ncp non-centrality parameter of the Beta distribution.
#' @param log,log.p logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.
#' @param lower.tail logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.
#' @return \samp{dbetagen} gives the density, \samp{pbetagen} gives the distribution function,
#' \samp{qbetagen} gives the quantile function, and \samp{rbetagen} generates random deviates.
#' @seealso \code{\link{Beta}}
#' @keywords distribution
#' @name betagen
#' @examples
#' curve(dbetagen(x, shape1=3, shape2=5, min=1, max=6), from = 0, to = 7)
#' curve(dbetagen(x, shape1=1, shape2=1, min=2, max=5), from = 0, to = 7, lty=2, add=TRUE)
#' curve(dbetagen(x, shape1=.5, shape2=.5, min=0, max=7), from = 0, to = 7, lty=3, add=TRUE)
#' @export
dbetagen <- function(x,shape1,shape2,min=0,max=1,ncp=0,log=FALSE)
{
  if(length(x) == 0) return(numeric(0))
  ow <- options(warn=-1)
  on.exit(options(ow), add = TRUE)
  min <- as.vector(min)
  max <- as.vector(max)
  x <- (x - min)/(max - min)
  if(missing(ncp))
  d <- dbeta(x, shape1=shape1, shape2=shape2) / (max-min)
  else
  {
  ncp<- as.vector(ncp)
  d <- dbeta(x, shape1=shape1, shape2=shape2, ncp=ncp) / (max-min)
  }
  if(log) d <- log(d)
  d[max <= min] <- NaN
  if(any(is.na(d))) warning("NaN in dbetagen")
  return(d)}

#' @rdname betagen
#' @export
pbetagen <- function(q,shape1,shape2,min=0,max=1,ncp=0,lower.tail = TRUE, log.p = FALSE)
{
  if(length(q) == 0) return(numeric(0))
  min <- as.vector(min)
  max <- as.vector(max)
  q2 <- (q - min)/(max-min)
  ow <- options(warn=-1)
  on.exit(options(ow), add = TRUE)
  if(missing(ncp))
  p <- pbeta(q2,shape1=shape1,shape2=shape2, lower.tail=lower.tail,log.p=log.p)
  else
  {
  ncp <- as.vector(ncp)
  p <- pbeta(q2,shape1=shape1,shape2=shape2, ncp=ncp, lower.tail=lower.tail, log.p=log.p)
  }
  # If min = max = q -> should return 1
  quel <- (abs(q - min) < (.Machine$double.eps^0.5)) &
          (abs(q - max) < (.Machine$double.eps^0.5))  #if min == max == q
  p[quel] <- if(lower.tail) 1 else 0
  if(log.p) p[quel] <- log(p[quel])

  p[max < min] <- NaN
  if(any(is.na(p))) warning("NaN in pbetagen")
  return(p)}

#' @rdname betagen
#' @export
qbetagen <- function(p,shape1,shape2,min=0,max=1,ncp=0,lower.tail=TRUE,log.p=FALSE)
{
  if(length(p) == 0) return(numeric(0))
  lout <- max(length(p),length(shape1),length(shape2),length(max),length(min),length(ncp))
  min <- rep(as.vector(min),length.out=lout)
  max <- rep(as.vector(max),length.out=lout)
  ow <- options(warn=-1)
  on.exit(options(ow), add = TRUE)
  if(missing(ncp))
    q <- qbeta(p,shape1=shape1, shape2=shape2, lower.tail=lower.tail, log.p=log.p)
  else
  {
    ncp <- as.vector(ncp)
    q <- qbeta(p,shape1=shape1, shape2=shape2, ncp=ncp, lower.tail=lower.tail, log.p=log.p)
  }
  q <- q * (max-min) + min
  q[max < min] <- NaN
  if(any(is.na(q))) warning("NaN in qbetagen")
  return(q)}

#' @rdname betagen
#' @export
rbetagen <- function(n,shape1,shape2,min=0,max=1,ncp=0)
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rbetagen")
  min <- rep(as.vector(min),length.out=n)
  max <- rep(as.vector(max),length.out=n)
  ow <- options(warn=-1)
  on.exit(options(ow), add = TRUE)
  if(missing(ncp))
    r <- rbeta(n, shape1=shape1, shape2=shape2)
  else
  {
    ncp <- as.vector(ncp)
    r <- rbeta(n, shape1=shape1, shape2=shape2, ncp=ncp)
  }
  r <- r*(max-min) + min
  r[max < min] <- NaN
  if(any(is.na(r))) warning("NaN in rbetagen")
  return(r)
}
