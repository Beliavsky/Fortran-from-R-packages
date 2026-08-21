##########################################################################################
## 1. Bivariate Poisson (BP)
##########################################################################################


#' @import Rcpp

dbp <- function(x, y, m0, m1, m2, log = FALSE, max = 500) {
  # max = 500: when counts exceed the max, p = 0
  if (m1 * m2 == 0) {  # previous code doesn't consider when mu1, mu2, or both is 0
    if ((x - y) %*% (m1 - m2) >= 0) {
      m <- min(x,y); s <- x + y - 2*m; mm <- max(m1, m2)
      result <- dpois (m, m0) * dpois (s, mm)
    } else { result <- 0 }
    if (log) {result <- log(result)}
    return(result)
  }
  f1 <- dpois(x, m1, log = TRUE)
  f2 <- dpois(y, m2, log = TRUE)
  m <- min(x,y); p <- m0/m1/m2; if (m == 0) {p <- 0} # in order not to make p NaN
  fun.a <- function(x, y, s, p, adj) {
    ifelse (p == 0, ifelse(s==0, 1, 0), exp(lchoose(x,s)+lchoose(y,s)+lfactorial(s)+ s*log(p) - adj))
  }
  if (max(x,y) > 100) {adj = 300} else {adj = 0}  # Handle numerical error for large numbers
  if (max(x,y) > max) {result <- ifelse(log, -Inf, 0)} else {
    f3 <- log(sum(sapply(0:m, function(s) fun.a(x=x, y=y, s=s, p=p, adj=adj)))) + adj
    result <- f1 + f2 - m0 + f3
    if (!log) {result <- exp(result)}
  }
  if (!is.finite(ifelse(log,exp(result),result))) {result <- 0}
  return(result)
}
dbp.vec <- Vectorize(dbp)

#' @rdname bp
#' @name bp
#' @aliases rbp
#' @aliases bp
#' @aliases lik.bp
#' 
#' @title The bivariate poisson distribution
#' 
#' @description random generation (\code{rbp}), maximum likelihood estimation (\code{bp}), 
#'    and log-likelihood. (\code{lik.bp})  for the bivariate Poisson 
#'    distribution with parameters equal to \code{(m0, m1, m2)}.
#' 
#' @param xvec,yvec a pair of bp random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(m0, m1, m2)}). 
#'    Either \code{param} or individual parameters (\code{m0, m1, m2}) 
#'    need to be provided.
#' @param m0,m1,m2 mean parameters of the Poisson variables. They must be positive.
#' @param n number of observations.
#' @param tol tolerance for judging convergence. \code{tol = 1e-8} by default.
#' 
#' @return 
#'  \itemize{
#'    \item \code{rbp} gives a pair of random vectors following BP distribution.
#'    \item \code{bp} gives the maximum likelihood estimates of a BP pair.
#'    \item \code{lik.bp} gives the log-likelihood of a set of parameters for a BP pair.
#'
#'  }
#'  
#' 
#' @examples
#' # generating a pair of random vectors
#' set.seed(1)
#' data1 <- rbp(n = 20, m0 = 1, m1 = 1, m2 = 1)
#' 
#' lik.bp(xvec = data1[, 1], yvec = data1[ ,2], 
#'           m0 = 1, m1 = 1, m2 = 1) 
#' 
#' bp(xvec = data1[,1], yvec = data1[,2])
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Liu, C., Preisser, J., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#'  Kocherlakota, S. & Kocherlakota, K. (1992). Bivariate Discrete Distributions. New York: Marcel Dekker.
#'  
#' @import Rcpp
#' @export
lik.bp <- function(xvec, yvec, m0, m1, m2, param=NULL) {
  if (!is.null(param)) {
    if (length(param) != 3) stop("length(param) must be 3.")
    m0 = param[1]; m1 = param[2]; m2 = param[3]
  }
  .check.m(c(m0, m1, m2))
  sum(log(do.call(dbp.vec, list(x = xvec, y = yvec, m0 = m0, m1 = m1, m2 = m2))))
}

#' @export
#' @rdname bp
rbp <- function(n, m0, m1, m2, param=NULL) {
  if (!is.null(param)) {
    if (length(param) != 3) stop("length(param) must be 3.")
    m0 = param[1]; m1 = param[2]; m2 = param[3]
  }
  if (length(n) != 1) {stop("length(n) must be 1.")}
  .check.m(c(m0, m1, m2))
  
  rmat <- matrix(rpois(n*3, lambda = c(m0, m1, m2)), n, 3, byrow=TRUE)
  xy <- rmat
  xy[,3] <- rmat[,1] + rmat[,3]
  xy[,2] <- rmat[,1] + rmat[,2]
  xy <- xy[,2:3]
  colnames(xy) <- c("x", "y")
  return(xy)
}

#' @export
#' @rdname bp
bp <- function(xvec, yvec, tol = 1e-6) {
  .check.input(xvec, yvec)
  
  xvec = as.integer(round(xvec, digits = 0))
  yvec = as.integer(round(yvec, digits = 0))
  
  len <- length(xvec)
  # xbar <- mean(xvec)
  # ybar <- mean(yvec)
  vec <- cbind(xvec,yvec)
  bar <- apply(vec,2,mean)
  m <- min(bar)
  lik.bp2 <- function(a) (-sum(apply(vec, 1, function(s) dbp(s[1], s[2], m0 = a, m1 = bar[1] - a, m2 = bar[2] - a, log = TRUE))))
  # result <- nlminb(m, lik, lower = 0)
  if (m == 0) {  # when either is all zero, then mu0 is automatically 0.
    result = c(par = 0, value = NA)  # lik not actually NA but can be obtained if necessary
  } else {
    result <- optim(m, lik.bp2, lower = 0, upper = m, method="Brent")
  }
  result <- c(mu0 = result$par, mu1 = bar[1] - result$par, mu2 = bar[2] - result$par, likelihood = -result$value)
  return(result)
}
if (FALSE) { # example
  bp(x,y)
  bp(x,z)
  set.seed(1000); a1 <- rpois(20,1); a2 <- rpois(20,2); a3 <- rpois(20,3)
  bp(a1+a2, a1+a3)
}


