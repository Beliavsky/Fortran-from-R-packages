##########################################################################################
## 2. Bivariate ZIP (BZIP)
##########################################################################################

## 2.1 basic functions for BvZIP
#' @import Rcpp

dbzip.a <- function(x, y = NULL, m0, m1, m2, p, log = FALSE) {
  fxy <- (1 - p) * dbp (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) + ifelse((x == 0 & y == 0), p,0)
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dbzip.a.vec <- Vectorize(dbzip.a)

#' @rdname bzip.a
#' @name bzip.a
#' @aliases rbzip.a
#' @aliases bzip.a
#' @aliases lik.bzip.a
#' 
#' @title The bivariate zero-inflated Poisson distribution (A)
#' 
#' @description random generation (\code{rbzip.a}), maximum likelihood estimation (\code{bzip.a}), 
#'    and log-likelihood. (\code{lik.bzip.a})  for the bivariate zero-inflated Poisson (A)
#'    distribution with parameters equal to \code{(m0, m1, m2, p)}.
#' 
#' @param xvec,yvec a pair of BZIP (A) random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(m0, m1, m2, p)}). 
#'    Either \code{param} or individual parameters (\code{m0, m1, m2, p}) 
#'    need to be provided.
#' @param m0,m1,m2  mean parameters of the Poisson variables. must be positive.
#' @param p zero-inflation probability 
#' @param n number of observations.
#' @param initial starting value of param for EM algorithm, a vector of nine values.
#' @param tol tolerance for judging convergence. \code{tol = 1e-8} by default.
#' @param showFlag if \code{TRUE}, the updated parameter estimates for each iteration 
#'   are printed out. If a positive integer, the updated parameter estimates for 
#'   iterations greater than \code{showFlag} are printed out.
#' 
#' @return 
#'  \itemize{
#'    \item \code{rbzip.a} gives a pair of random vectors following BZIP (A) distribution.
#'    \item \code{bzip.a} gives the maximum likelihood estimates of a BZIP (A) pair.
#'    \item \code{lik.bzip.a} gives the log-likelihood of a set of parameters for a BZIP (A) pair.
#'
#'  }
#'  
#' 
#' @examples
#' # generating a pair of random vectors
#' set.seed(1)
#' data1 <- rbzip.a(n = 20, m0 = 1, m1 = 1, m2 = 1, p = 0.5)
#' 
#' lik.bzip.a(xvec = data1[, 1], yvec = data1[ ,2], 
#'           m0 = 1, m1 = 1, m2 = 1, p = 0.5)
#' 
#' bzip.a(xvec = data1[,1], yvec = data1[,2], showFlag = FALSE)
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Liu, C., Preisser, J., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#'  Li, C. S., Lu, J. C., Park, J., Kim, K., Brinkley, P. A., & Peterson, J. P. (1999). 
#'  Multivariate zero-inflated Poisson models and their applications. Technometrics, 41, 29-38.
#' 
#' @import Rcpp
#' @export
#' 
lik.bzip.a <- function(xvec, yvec, m0, m1, m2, p, param=NULL) {
  if (!is.null(param)) {
    if (length(param) != 4) stop("length(param) must be 4.")
    m0 = param[1]; m1 = param[2]; m2 = param[3]; p = param[4]
  }
  .check.input(xvec, yvec)
  .check.m(c(m0, m1, m2))
  .check.p(c(p, 1 - p, 0, 0))
  sum(dbzip.a.vec(x = xvec, y = yvec, m0 = m0, m1 = m1, m2 = m2, p = p, log=TRUE))
}

#' @export
#' @rdname bzip.a
rbzip.a <- function(n, m0, m1, m2, p, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 4) stop("length(param) must be 4.")
    m0 = param[1]; m1 = param[2]; m2 = param[3]; p = param[4]
  }
  if (length(n) != 1) {stop("length(n) must be 1.")}
  .check.m(c(m0, m1, m2))
  .check.p(c(p, 1 - p, 0, 0))
  
  rmat <- matrix(rpois(n*3, lambda = c(m0, m1, m2)), n, 3, byrow=TRUE)
  xy <- rmat
  xy[,3] <- rmat[,1] + rmat[,3]
  xy[,2] <- rmat[,1] + rmat[,2]
  xy <- xy[,2:3]
  colnames(xy) <- c("x", "y")
  z <- rbinom(n, 1, p)
  xy <- xy * (1 - z)
  return(xy)
}

## 2.2 param estimation
# formal EM algorithm
#' @export
#' @rdname bzip.a
#' @importFrom stats dpois
bzip.a <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE) { #MLE based on score equations : fail (not convex)
  .check.input(xvec, yvec)
  
  if (!is.null(initial)) {
    if (length(initial) != 4) {stop("length(initial) must be 4.")}
    .check.m(initial[1:3])
    .check.p(c(initial[4], 1 - initial[4], 0, 0))
  }
  
  # counter, likelihood, param.prev for recursive function
  len <- length(xvec)          # n
  vec <- data.frame(xvec = xvec, yvec = yvec)
  sum.x.y <- apply(vec, 2, sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(c(mu0 = NA, mu1 = NA, mu2 = NA, p= NA))} # everything is zero ==> nonestimable, set pi = 0
  
  # E-step
  fun.cond.exp <- function(x, y, m0, m1, m2, p) {
    if (x+y == 0) {
      cond.prob <- p / (p + (1-p) * exp(-m0 -m1 -m2))
      cond.expt <- c(1, m0, m1, m2) * cond.prob
    } else {
      m = min(x,y)
      prob <- sapply(0:m, function(u) {dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)})
      if (sum(prob) == 0) {
        # using the ratio of two probs
        if (m==0) {prob <- 1} else {
          prob <- sapply(0:(m-1), function(u) (m0/m1/m2*(x-u)*(y-u)/(u+1)))
          prob <- cumprod(c(1,prob))
          #prob <- sapply(0:m, function(u) {dbp(x = u, y=0, m0 = 0, m1 = m0, m2 = 0) * dbp(x = x-u, y=0, m0 = 0, m1 = m1, m2 = 0) *
          #    dbp(x = 0, y= y-u, m0 = 0, m1 = 0, m2 = m2)*10^100})
          # used dbp instead of simply dpois to handle large numbers (when mu=1, x = 172, dpois=0 not small number)
          # and adjust by 10^100 (num/den cancel out adj.)
        }
      }
      eu <- sum((0:m)*prob)/sum(prob)
      cond.expt <- c(0, 0, x, y) + eu *c(0, 1, -1, -1)
    }
    return(cond.expt[c(2:4,1)]) # changing (p, m0, m1, m2) to (m0, m1, m2, p) 
  }
  fun.cond.exp <- Vectorize(fun.cond.exp)
  
  # M-step
  param.update <- function(x, y, m0, m1, m2, p) {
    result <- fun.cond.exp(x = x, y = y, m0 = m0, m1 = m1, m2 = m2, p = p)
    return(apply(result, 1, mean))
  }
  
  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- c(1, 1, 1, 0.5)
  }
  
  # Repeat
  iter = 0
  param = initial
  repeat {
    iter = iter + 1
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec, m0 = param[1], m1 = param[2], m2 = param[3], p = param[4])
    if (showFlag) {message("iter ", iter, ", param ", 
                           paste0(c("mu0 ", "mu1 ", "mu2 ", "p "), round(param, 5), ", "), 
                           "lik ", round(lik.bzip.a(xvec, yvec, param = param), 4))}
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("mu0", "mu1", "mu2", "p")
      return(param)
      break
    }
  }
  return(param)
}

moment.bzip <- function(m0, m1, m2, p) {
  mean.bp <- (m0 + c(m1, m2))
  var.bp <- diag(mean.bp); var.bp[c(2,3)] <- m0
  mean <- (1 - p)*mean.bp
  var <- (1 - p)*var.bp + p*(1 - p) * mean.bp %o% mean.bp
  cor <- var[2] / sqrt(prod(diag(var)))
  return(list(mean = mean, var = var, cor = cor))
}


##########################################################################################
## 3. Bivariate zip.b: General bzip with marginal ZIP condition (6params)
##########################################################################################
## 3.1 basic functions for bzip.b


#' @import Rcpp

dbzip.b <- function(x, y = NULL, m0, m1, m2, p1, p2, p3, p4 = 1 - p1 - p2 - p3, log = FALSE) {
  fxy <- p1 * dbp (x=x, y=y, m0 = m0, m1 = m1, m2 = m2) +
    p2 * {if (y == 0) dbp (x=x, y=y, m0 = 0, m1 = m0 + m1, m2 = 0) else 0} +
    p3 * {if (x == 0) dbp (x=x, y=y, m0 = 0, m1 = 0, m2 = m0 + m2) else 0} +
    p4 * {if (x + y == 0) 1 else 0}
  if (log) {fxy <- log(fxy)}
  return(fxy)
}
dbzip.b.vec <- Vectorize(dbzip.b)

#' @rdname bzip.b
#' @name bzip.b
#' @aliases rbzip.b
#' @aliases bzip.b
#' @aliases lik.bzip.b
#' 
#' @title The bivariate zero-inflated Poisson distribution (B)
#' 
#' @description random generation (\code{rbzip.b}), maximum likelihood estimation (\code{bzip.b}), 
#'    and log-likelihood. (\code{lik.bzip.b})  for the bivariate zero-inflated Poisson (B)
#'    distribution with parameters equal to \code{(m0, m1, m2, p1, p2, p3, p4)}.
#' 
#' @param xvec,yvec a pair of BZIP (B) random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(m0, m1, m2, p1, p2, p3, p4)}). 
#'    Either \code{param} or individual parameters (\code{m0, m1, m2, p1, p2, p3, p4}) 
#'    need to be provided.
#' @param m0,m1,m2  mean parameters of the Poisson variables. must be positive.
#' @param p1,p2,p3,p4 proportions summing up to 1 (\code{p1 + p2 + p3 + p4 = 1}). 
#' \code{p1} is the probability of both latent Poisson variables being observed. 
#' \code{p2} is the probability of only the first Poisson variables being observed.
#' \code{p3} is the probability of only the second Poisson variables being observed, and
#' \code{p4} is the probability of both Poisson variables being dropped out. 
#' @param n number of observations.
#' @param initial starting value of param for EM algorithm, a vector of nine values.
#' @param tol tolerance for judging convergence. \code{tol = 1e-8} by default.
#' @param maxiter maximum number of iterations allowed. \code{tol = 50000} by default.
#' @param showFlag if \code{TRUE}, the updated parameter estimates for each iteration 
#'   are printed out. If a positive integer, the updated parameter estimates for 
#'   iterations greater than \code{showFlag} are printed out.
#' 
#' @return 
#'  \itemize{
#'    \item \code{rbzip.b} gives a pair of random vectors following BZIP (B) distribution.
#'    \item \code{bzip.b} gives the maximum likelihood estimates of a BZIP (B) pair.
#'    \item \code{lik.bzip.b} gives the log-likelihood of a set of parameters for a BZIP (B) pair.
#'
#'  }
#'  
#' 
#' @examples
#' # generating a pair of random vectors
#' set.seed(1)
#' data1 <- rbzip.b(n = 20, m0 = 1, m1 = 1, m2 = 1, 
#'                 p1 = 0.5, p2 = 0.2, p3 = 0.2, p4 = 0.1)
#' 
#' lik.bzip.b(xvec = data1[, 1], yvec = data1[ ,2], 
#'           m0 = 1, m1 = 1, m2 = 1, 
#'           p1 = 0.5, p2 = 0.2, p3 = 0.2, p4 = 0.1)
#' 
#' bzip.b(xvec = data1[,1], yvec = data1[,2], showFlag = FALSE)
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Liu, C., Preisser, J., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#' @import Rcpp
#' @export
#' 
lik.bzip.b <- function(xvec, yvec, m0 , m1, m2, p1, p2, p3, p4, param = NULL){
  if (!is.null(param)) {
    if (length(param) != 7) stop("length(param) must be 7.")
    m0 = param[1]; m1 = param[2]; m2 = param[3]
    p1 = param[4]; p2 = param[5]; p3 = param[6]; p4 = param[7]
  }
  sum(dbzip.b.vec(x = xvec, y = yvec, m0 = m0, m1 = m1, m2 = m2, 
                  p1 = p1, p2 = p2, p3 = p3, p4 = p4, log = TRUE))
}

#' @export
#' @rdname bzip.b
rbzip.b <- function(n, m0, m1, m2, p1, p2, p3, p4, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 7) stop("length(param) must be 4.")
    m0 = param[1]; m1 = param[2]; m2 = param[3]; p1 = param[4]
    p2 = param[5]; p3 = param[6]; p4 = param[7]
  }
  if (length(n) != 1) {stop("length(n) must be 1.")}
  .check.m(c(m0, m1, m2))
  .check.p(c(p1, p2, p3, p4))
  
  rmat <- matrix(rpois(n*3, lambda = c(m0, m1, m2)), n, 3, byrow=TRUE)
  xy <- rmat
  xy[,3] <- rmat[,1] + rmat[,3]
  xy[,2] <- rmat[,1] + rmat[,2]
  xy <- xy[,2:3]
  colnames(xy) <- c("x", "y")
  E <- t(rmultinom(n, 1, c(p1, p2, p3, p4)))
  z <- cbind(E[,1] + E[,2], E[,1] + E[,3])
  
  xy <- xy * z
  return(xy) 
}

## 3.2 parameter estimation
# Method of moment, for starting values of MLE
mme.bzip.b <- function(xvec, yvec) {
  m.x <- mean(xvec); m.x2 <- mean(xvec^2)
  m.y <- mean(yvec); m.y2 <- mean(yvec^2)
  m.x2.y <- mean(xvec^2*yvec); m.y2.x <- mean(yvec^2*xvec)
  m.x.y <- mean(xvec*yvec)
  w <- (m.x2.y + m.y2.x) / (2 * m.x.y)
  lambda1 <- m.x2 / m.x - 1
  lambda2 <- m.y2 / m.y - 1
  mu0 <- (lambda1*lambda2) * (w - (2 + lambda1 + lambda2)/2) / (lambda1 + lambda2 + 1 - w)
  mu0 <- max(1e-10, min(mu0, m.x, m.y))
  mu1 <- m.x - mu0
  mu2 <- m.y - mu0
  
  pi1 <- m.x.y / (lambda1*lambda2 + mu0)
  pi2 <- m.x / lambda1 - pi1
  pi3 <- m.y / lambda2 - pi1
  pi4 <- 1 - pi1 - pi2 - pi3
  pi <- c(pi1, pi2, pi3, pi4)
  pi <- pmax(0, pmin(pi,1))
  pi <- pi/sum(pi)
  return(data.frame(mu0 = mu0, mu1 = mu1, mu2 = mu2,
                    pi1 = pi[1], pi2 = pi[2], pi3 = pi[3], pi4 = pi[4]))
}

#' @export
#' @rdname bzip.b
#' @importFrom stats dpois
bzip.b <- function(xvec, yvec, tol = 1e-6, initial = NULL, showFlag = FALSE, maxiter = 200) {
  .check.input(xvec, yvec)
  if (!is.null(initial)) {
    if (length(initial) != 7) {stop("length(initial) must be 7.")}
    .check.m(initial[1:3])
    .check.p(initial[6:9])
  }
  # counter, likelihood, param.prev for recursive function
  # initial guess manipulation (for small counts pi4=0)
  len <- length(xvec)          # n
  vec <- data.frame(xvec=xvec, yvec=yvec)
  sum.x.y <- apply(vec,2,sum)  # mu0+mu1, mu0+mu2
  if (sum(sum.x.y) ==0 ) { return(c(mu0 = NA, mu1 = NA, mu2 = NA, p1 = NA, 
                                             p2 = NA, p3 = NA, p4 = NA))} # everything is zero ==> nonestimable
  
  # E-step
  fun.cond.exp.a <- function(x, y, m0, m1, m2, p1, p2, p3, p4) {
    pra <- p1 * dbp(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    prb <- p2 * (y==0) * dbp(x=x, y=0, m0 =0, m1 = m0 + m1, m2 = 0)
    prc <- p3 * (x==0) * dbp(x=0, y=y, m0 =0, m1 = 0, m2 = m0 + m2)
    prd <- p4 * (x + y ==0)
    pragr <- c(pra, prb, prc, prd)
    prsum <- sum(pragr)
    
    ### Conditional expectations (EE) for all profile cases of (x, y)
    EE <- (pragr/prsum)
    
    ### Conditional expectations (EkU) for all profile cases of (x, y)
    m = min(x,y)
    num <- sum(sapply(0:m, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1) * dpois(x = y - u, lambda = m2)}))
    den <- dbp(x=x, y=y, m0 = m0, m1 = m1, m2 = m2)
    EE1U <- EE[1] * num/den
    
    num <- sum(sapply(0:x, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = x - u, lambda = m1)}))
    den <- dbp(x=x, y=0, m0 = 0, m1 = m0 + m1, m2 = 0)
    EE2U <- EE[2] * num/den
    
    num <- sum(sapply(0:y, function(u) {u * dpois(x = u, lambda = m0) * dpois(x = y - u, lambda = m2)}))
    den <- dbp(x=0, y=y, m0 = 0, m1 = 0, m2 = m0 + m2)
    EE3U <- EE[3] * num/den
    
    EE4U <- EE[4] * m0
    
    m0.den <- 1
    m0.num <- (EE1U + EE2U + EE3U + EE4U)
    m1.den <- sum(EE[1:2])
    m1.num <- m1.den * x - EE1U - EE2U
    m2.den <- sum(EE[c(1,3)])
    m2.num <- m2.den * y - EE1U - EE3U
    
    return(c(EE, m0.den, m0.num, m1.den, m1.num, m2.den, m2.num))
  }
  fun.cond.exp <- Vectorize(fun.cond.exp.a)
  
  # M-step
  param.update <- function(x, y, m0, m1, m2, p1, p2, p3, p4) {
    result <- fun.cond.exp(x, y, m0, m1, m2, p1, p2, p3, p4)
    result[result < 0] <- 0
    result <- apply(result,1,mean)
    result2 <- c(0,0,0, result[1:4])
    result2[1] <- result[6]/result[5]
    result2[2] <- result[8]/result[7]
    result2[3] <- result[10]/result[9]
    result2[is.na(result2)] <- as.numeric(c(m0, m1, m2, p1, p2, p3, p4))[is.na(result2)]
    return(result2)
  }
  
  # initial guess
  if (is.null(initial)) { # when initial(starting point) is not provided
    initial <- mme.bzip.b(xvec = xvec, yvec = yvec)
    if (sum(is.na(initial)) > 0) {
      initial <- rep(0, 7)
      initial[4:7] <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
      initial <- initial/sum(initial)      # relative freq
      initial[1:3] <- c(0.001, sum.x.y/len/(1-initial[7]))
    }
    if (min(sum.x.y) < 5 & min(sum.x.y) > 0) {
      initial[5:6] <- (initial[7] - 1e-10) * sum.x.y/sum(sum.x.y)
      initial[7] <- 1e-10}
  }
  initial <- pmax(initial, rep(1e-10, 7))
  
  # Repeat
  iter = 0
  param = initial
  repeat {
    if (iter >= maxiter) { warning("EM exceeded maximum number of iterations")
      param <- rep(NA,7)
      names(param) <- c("mu0", "mu1", "mu2", "p1", "p2", "p3", "p4")
      return(param)}
    
    iter = iter + 1
    param.old <- param # saving old parameters
    param <- param.update (x = xvec, y = yvec,m0=param[1], m1=param[2], m2=param[3],
                           p1 = param[4], p2 = param[5], p3 = param[6], p4 = param[7])
    if (showFlag) {message("iter ", iter, ", param ", 
                           paste0(c("mu0 ", "mu1 ", "mu2 ", "p1 ", "p2 ", "p3 ", "p4 "), round(param, 5), ", "), 
                           "lik ", round(lik.bzip.b(xvec, yvec, param = param), 4))}
    if (max(abs(param - param.old)) <= tol) {
      names(param) <- c("mu0", "mu1", "mu2", "p1", "p2", "p3", "p4")
      return(param)
      break
    }
  }
  return(param)
}