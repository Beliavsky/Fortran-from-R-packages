##########################################################################################
## 4. Bivariate Zero-Inflated Negative Binomial (BZINB)
##########################################################################################


#' @import Rcpp stats
###@importFrom stats rpois rgamma rmultinom rbinom var cor optim
###@importFrom stats dpois dgamma dmultinom dbinom qlogis
###@importFrom stats cor dnbinom dpois optim qlogis rbinom rgamma rmultinom rpois setNames var

dbzinb <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, log=FALSE) {
  dxy <- dbnb(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a2, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(ifelse(log, log(result), result))
}
dbzinb.vec <- Vectorize(dbzinb)

# weight of nondropout
wbzinb.i <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 9) stop("length(param) must be 9.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
    p1 = param[6]; p2 = param[7]; p3 = param[8]; p4 = param[9]
  }
  
  dxy <- dbnb(x=x, y=y, a0=a0, a1=a1, a2=a2, b1=b1, b2=b2, log=FALSE)
  dx <- dnbinom(x=x, a0+a1, 1/(1+b1))
  dy <- dnbinom(x=y, a0+a2, 1/(1+b2))
  result <- dxy * p1 + dx * ifelse(y==0,p2,0) + dy * ifelse(x==0,p3,0) + ifelse(x+y==0,p4,0)
  return(dxy * p1 / result)
}
wbzinb <- Vectorize(wbzinb.i, vectorize.args = c("x", "y"))

weighted.cor.base <- function(xvec, yvec, weight) {
  Exy = weighted.mean(xvec * yvec, w = weight)
  Exx = weighted.mean(xvec * xvec, w = weight)
  Eyy = weighted.mean(yvec * yvec, w = weight)
  Ex = weighted.mean(xvec, w = weight)
  Ey = weighted.mean(yvec, w = weight)
  cor = (Exy - Ex * Ey) / sqrt((Exx - Ex^2) * (Eyy - Ey^2))
  return(cor)
}

#' @rdname weighted.pc
#' @name weighted.pc
#' 
#' @title Weighted Pearson Correlation (WPC) based on bivariate zero-inflated 
#' negative binomial (BZINB) model
#' 
#' @description weighted.pc calculates Pearson's correlation with less weights for pairs 
#' containing zero(s). The weights are determined by BZINB model.
#'  
#' @param xvec,yvec a pair of bzinb random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(a0, a1, a2, b1, b2, p1, p2, p3, p4)}). 
#'    See \code{bzinb} for details.
#'    If \code{param} is \code{null}, it will be estimated by \code{bzinb()}.
#' @param ... optional arguments used passed to \code{bzinb}, when \code{param} is \code{null}.
#' 
#' @return weighted Pearson correlation (WPC) estimate
#' 
#' @examples
#' # generating a pair of random vectors
#' set.seed(2)
#' data1 <- rbzinb(n = 20, a0 = 1, a1 = 1, a2 = 1, 
#'                 b1 = 1, b2 = 1, p1 = 0.5, p2 = 0.2, 
#'                 p3 = 0.2, p4 = 0.1)
#' 
#' weighted.pc(xvec = data1[,1], yvec = data1[,2], 
#'             param = c(0.769, 0.041, 0.075, 3.225, 1.902, 0.5, 0.084, 1e-20, 0.416))
#' weighted.pc(xvec = data1[,1], yvec = data1[,2])
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Preisser, J., Liu, C., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#'  
#' @import Rcpp
#' @export
#' @useDynLib bzinb
#' 
weighted.pc <- function(xvec, yvec, param = NULL, ...) {
  .check.input(xvec, yvec)
  
  if (!is.null(param)) {
    if (length(param) != 9) stop("length(param) must be 9.")
    .check.ab(param[1:5])
    .check.p(param[6:9])
  } else {
    bzinb.out = bzinb(xvec = xvec, yvec = yvec) 
    param = bzinb.out$coefficients[, 1]
  }
  
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  
  weight = wbzinb(xy.reduced$x, xy.reduced$y, param = param)
  wcor = weighted.cor.base(xy.reduced$x, xy.reduced$y, weight = xy.reduced$freq * weight)
  wcor
}

#' @rdname bzinb
#' @name bzinb
#' @aliases rbzinb
#' @aliases bzinb
#' @aliases lik.bzinb
#' 
#' @title The bivariate zero-inflated negative binomial distribution
#' 
#' @description random generation (\code{rbzinb}), maximum likelihood estimation (\code{bzinb}), 
#'    and log-likelihood. (\code{lik.bzinb})  for the bivariate zero-inflated negative binomial 
#'    distribution with parameters equal to \code{(a0, a1, a2, b1, b2, p1, p2, p3, p4)}.
#' 
#' @param xvec,yvec a pair of bzinb random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(a0, a1, a2, b1, b2, p1, p2, p3, p4)}). 
#'    Either \code{param} or individual parameters (\code{a0, a1, a2, b1, b2, p1, p2, p3, p4}) 
#'    need to be provided.
#' @param a0,a1,a2 shape parameters of the latent gamma variables. They must be positive.
#' @param b1,b2 scale parameters for the latent gamma variables. They must be positive.
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
#' @param vcov if \code{TRUE}, the variance-covariance matrix and information matrix 
#'   are returned.
#' 
#' @return 
#'  \itemize{
#'    \item \code{rbzinb} gives a pair of random vectors following BZINB distribution.
#'    \item \code{bzinb} gives the maximum likelihood estimates of a BZINB pair.
#'      \itemize{
#'        \item \code{rho} estimate and standard error of the underlying correlation (\eqn{\rho}) and (\eqn{logit(\rho)})
#'        \item \code{coefficients} estimate and standard error of the BZINB parameters
#'        \item \code{lik} log-likelihood of the maximum likelihood estimate
#'        \item \code{iter} total number of iterations
#'        \item \code{info} information matrix. Provided when \code{vcov} is \code{TRUE}.
#'        \item \code{vcov} variance-covariance matrix. Provided when \code{vcov} is \code{TRUE}.
#'      }
#'    \item \code{lik.bzinb} gives the log-likelihood of a set of parameters for a BZINB pair.
#'
#'  }
#'  
#' @details
#' EM theoretically guarantees higher likelihood at each iteration than that of previous iterations. 
#' See Dempster, Laird, and Rubin (1977). This guarantee comes with an assumption that there is no numerical
#' error in conditional likelihood maximization at each iteration. Small errors can cause decreasing likelihood 
#' especially when the iterations reach the point of convergence. Due to this technical error, the EM continues after it reaches
#' the maximum likelihood point (up to 100 iterations). However, the final estimate being returned is the parameter values at
#' the maximum likelihood.
#'  
#' 
#' @examples
#' # generating a pair of random vectors
#' set.seed(2)
#' data1 <- rbzinb(n = 100, a0 = 2, a1 = 1, a2 = 1, 
#'                 b1 = 1, b2 = 1, p1 = 0.5, p2 = 0.2, 
#'                 p3 = 0.2, p4 = 0.1)
#' 
#' lik.bzinb(xvec = data1[, 1], yvec = data1[ ,2], 
#'           a0 = 1, a1 = 1, a2 = 1, b1 = 1, b2 = 1, 
#'           p1 = 0.5, p2 = 0.2, p3 = 0.2, p4 = 0.1)
#' 
#' bzinb(xvec = data1[,1], yvec = data1[,2], showFlag = FALSE)
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Preisser, J., Liu, C., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#'  Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from 
#'  incomplete data via the EM algorithm. Journal of the Royal Statistical Society: 
#'  Series B (Methodological), 39(1), 1-22.
#'  
#' @import Rcpp
#' @export
#' @useDynLib bzinb
#' 
lik.bzinb <- function(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, param=NULL) {
  if (!is.null(param)) {
    if (length(param) != 9) stop("length(param) must be 9.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
    p1 = param[6]; p2 = param[7]; p3 = param[8]; p4 = param[9]
  }
  .check.input(xvec, yvec)
  .check.ab(c(a0, a1, a2, b1, b2))
  .check.p(c(p1, p2, p3, p4))
  sum(log(do.call(dbzinb.vec, list(x = xvec, y = yvec, a0 = a0, a1 = a1, a2 = a2, 
                                   b1 = b1, b2 = b2, p1 = p1, p2 = p2, p3 = p3, p4 = p4))))
}


#' @export
#' @rdname bzinb
rbzinb <- function(n, a0, a1, a2, b1, b2, p1, p2, p3, p4, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 9) stop("length(param) must be 9.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
    p1 = param[6]; p2 = param[7]; p3 = param[8]; p4 = param[9]
  }
  if (length(n) != 1) {stop("length(n) must be 1.")}
  .check.ab(c(a0, a1, a2, b1, b2))
  .check.p(c(p1, p2, p3, p4))
  
  rmat <- matrix(rgamma(n*3, shape = c(a0, a1, a2), rate = 1/b1), n, 3, byrow=TRUE)
  rmat2 <- rmat
  rmat2[,3] <- rmat[,1] + rmat[,3]
  rmat2[,2] <- rmat[,1] + rmat[,2]
  rmat2 <- rmat2[,2:3]
  rmat2[,2] <- rmat2[,2]*b2/b1
  uv <- matrix(rpois(n*2, rmat2), n, 2)

  E <- t(rmultinom(n, 1, c(p1, p2, p3, p4)))
  z <- cbind(E[,1]+E[,2], E[,1]+E[,3])

  xy <- uv * z
  colnames(xy) <- c("x", "y")

  return(xy)
}

trueCor <- function(a0, a1, a2, b1, b2) a0 * sqrt(b1 * b2 /(b1+1) /(b2+1) /(a0 + a1) /(a0 + a2))


### 2.em

# # not used
# dbzinb.expt.vec <- function(xvec, yvec, freq, n = length(freq), a0, a1, a2, b1, b2, p1, p2, p3, p4) {
#   .check.input()
#   
#   result <- rep(0, 12)
#   dbzinb_expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
#   names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
#   result
# }
# 
# # not used
# # step-by-step expt calculation for information matrix
# dbzinb.expt.se <- function(x, y, a0, a1, a2, b1, b2, p1, p2, p3, p4) {
#   result <- rep(0, 12)
#   dbzinb_expt(x, y, freq = 1, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
#   names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
#   result
# }
# dbzinb.expt.se.mat <- Vectorize(dbzinb.expt.se)
# 
# dbzinb.expt.mat <- function(xvec, yvec, freq, n = sum(freq), a0, a1, a2, b1, b2, p1, p2, p3, p4) {
#   result <- rep(0, 12)
#   dbzinb_expt_vec(xvec, yvec, freq, n, a0, a1, a2, b1, b2, p1, p2, p3, p4, result)
#   names(result) <- c("logdensity", paste0("R", 0:2, ".E"), paste0("log.R", 0:2, ".E"), paste0("E",1:4,".E"), "v.E")
#   result
# }
# 
# if (FALSE) {
#   tmp <- dbzinb.expt.vec(c(1,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
#   tmp <- dbzinb.expt.vec(c(0,1,1),c(0,1,2),1,1,1,1,2,.25,.25,.25,.25)
#   tmp <- dbzinb.expt.vec(extractor(1),extractor(2),1,1,1,1,2,.25,.25,.25,.25)
#   t(tmp)[21:40,]
#   dbzinb.expt.vec(c(10,1,2),c(10,1,1), 1.193013282, 0.003336139, 0.002745513, 3.618842924, 3.341625901, .25,.25,.25,.25)
# }

abp.names <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4") # global variable
expt.names <- c("lik", "ER0", "ER1", "ER2", "ElogR0", "ElogR1", "ElogR2", "EE1", "EE2", "EE3", "EE4", "EV")

bzinb.base <- function (xvec, yvec, initial = NULL, tol = 1e-8, maxiter=50000, showFlag=FALSE, vcov = FALSE,
                        bnb = 0) {
  se = TRUE  # se is estimated by default
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n.reduced <- as.integer(length(xy.reduced$freq))
  

  if (max(xvec)==0 & max(yvec)==0) {return(c(rep(1e-10,5),1,0,0,0, 0, 1, 0, if (se) {rep(NA, 11)}))} # 9 params, lik, iter, pureCor, and 11 se's
  info <- if (se) {matrix(0, ncol = 8, nrow = 8, dimnames = list(abp.names[-9], abp.names[-9]))} else {0}

  # initial guess
  if (is.null(initial)) {
    xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
    s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)|is.na(s2.y)) {s2.x <- s2.y <- 1}
    cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
    zero <- sum(xvec == 0 & yvec == 0) / n

    initial <- rep(NA,9)
    names(initial) <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4")
    initial[4] <- s2.x /ifelse(xbar==0,1e-4, xbar)
    initial[5] <- s2.y /ifelse(ybar==0,1e-4, ybar)
    initial[2:3] <- c(xbar,ybar)/pmax(initial[4:5], c(0.1,0.1))
    initial[1] <- min(initial[2:3]) * abs(cor.xy)
    initial[2:3] <-  initial[2:3] - initial[1]

    initial[6:9] <- bin.profile(xvec, yvec)   # freq of each zero-nonzero profile
    initial[6:9] <- initial[6:9]/sum(initial[6:9])      # relative freq
    initial <- pmax(initial, 1e-5)
    if(is.na(sum(initial))) { initial[is.na(initial)] <- 1}
  } else {
    names(initial) <- abp.names
  }
  ## tmp.initial <<- initial
  param = initial
  lik = -Inf
  lik.vec = rep(0, maxiter + 1)
  nonconv = 0L

  em.out <- em(param2 = param, xvec = xy.reduced$x, yvec = xy.reduced$y, 
       freq = xy.reduced$freq, n = n.reduced, se = as.integer(se), 
       maxiter = as.integer(maxiter), tol = as.double(tol), 
       showFlag = as.integer(showFlag), bnb = bnb)
  names(em.out) <- c("param2",              # "xvec", "yvec", "freq", "n", 
                     "expt", "info",        # "se", 
                     "iter",                # "maxiter", "tol", "showFlag", 
                     "nonconv", "trajectory") # , "bnb")
  
  # overwriting the original object
  param = setNames(em.out$param2, abp.names)
  expt  = setNames(em.out$expt, expt.names)
  info = if (se) {matrix(em.out$info, 8, 8, 
                         dimnames = list(abp.names[-9], abp.names[-9]))} else 0
  iter  = em.out$iter
  nonconv = em.out$nonconv
  trajectory = em.out$trajectory
  
  if (nonconv == 1) warning("The iteration exited before reaching convergence.")
  
  # tmp.expt <<- expt
  # tmp.traj <<- lik.vec
  # underlying correlation (rho)
  rho <- param[1]/sqrt((param[1] + param[2]) * (param[1] + param[3])) *
    sqrt(param[4] * param[5] /(param[4] + 1) /(param[5] + 1))
  logit.rho <- qlogis(rho)
  
  if (bnb) {info = info[1:5, 1:5]}
  
  if (se) {
    par.names <- if (bnb) {abp.names[1:5]} else {abp.names}
    dim.std <- if (bnb) {7} else {11}
    fullrank <- if (bnb) {5} else {8}
    
    qr.info <- try(qr(info))
    if (class(qr.info)[1] == "try-error") {
      warning ("The information matrix has NA/NaN/Inf and thus the standard error is not properly estimatd.")
      std.param = setNames(rep(NA, dim.std), c(par.names, "rho", "logit.rho"))
      cov.mat <- NA
    } else if (qr(info)$rank < fullrank) {
      warning ("The information matrix is (essentially) not full rank, and thus the standard error is not reliable.")
      std.param = setNames(rep(NA, dim.std), c(par.names, "rho", "logit.rho"))
      cov.mat <- NA
    } else {
      cov.mat <- try(solve(info))
      if (class(cov.mat)[1] == "try-error") {
        std.param = setNames(rep(NA, dim.std), c(par.names, "rho", "logit.rho"))
        cov.mat <- NA
      } else {
        # variance of p4 hat
        if (!bnb) var.p4 <- sum (cov.mat[6:8, 6:8]) # = sum_i,j cov(pi, pj)
        
        # variance of rho hat
        # rho <- a0/sqrt((a0 + a1) * (a0 + a2)) *sqrt(b1 *b2 /(b1 + 1) /(b2 + 1))
        
        # d.g <- rho * c(1/a0 - 1/{2*(a0 + a1)} - 1/{2*(a0 + a2)}, - 1/{2*(a0 + a1)}, - 1/{2*(a0 + a2)},
        #                1/{2 *b1 *(b1 + 1)}, 1/{2 *b2 *(b2 + 1)})
        d.g <- rho * c(1/param[1] - 1/{2*(param[1] + param[2])} - 1/{2*(param[1] + param[3])}, 
                       - 1/{2*(param[1] + param[2])}, 
                       - 1/{2*(param[1] + param[3])},
                       1/{2 *param[4] *(param[4] + 1)}, 
                       1/{2 *param[5] *(param[5] + 1)})
        
        var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g
        
        # variance of logit(rho hat)
        var.logit.rho <- var.rho / rho^2 / (1-rho)^2
        # std.param = sqrt(c(setNames(diag(cov.mat), abp.names[1:8]), 
        #                    p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
        if (!bnb) {
          std.param = sqrt(c(diag(cov.mat), 
                 p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
        } else {
          std.param = sqrt(c(diag(cov.mat), rho=var.rho, logit.rho = var.logit.rho))
        }
        
      }
    } 
    
  }
 
  result <- list(rho = matrix(c(rho, logit.rho, if(se) std.param[c("rho", "logit.rho")] else rep(NA, 2)),
                              ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                 coefficients = matrix(c(param, if(se) std.param[1:9] else rep(NA, 9)),
                                       ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                 lik = expt[1],
                 iter = iter)
  if (se & vcov) {
    result$info = info
    result$vcov = cov.mat
  }
  return(result)
}

#' @export
#' @rdname bzinb
bzinb <- function(xvec, yvec, initial = NULL, tol = 1e-8, maxiter = 50000, showFlag = FALSE,
                  vcov = FALSE) {
  .check.input(xvec, yvec)
  if (!is.null(initial)) {
    if (length(initial) != 9) {stop("length(initial) must be 9.")}
    .check.ab(initial[1:5])
    .check.p(initial[6:9])
  }
  
  # if (!is.integer(xvec) | !is.integer(yvec)) stop("xvec and yvec must be integers.")
  # nonnegative
  # len(xvec) == len(yvec)
  # any(is.na(xvec))
  xvec = as.integer(round(xvec, digits = 0))
  yvec = as.integer(round(yvec, digits = 0))
  result <- try(bzinb.base(xvec, yvec, initial = initial, tol = tol, maxiter = maxiter, 
                           showFlag = showFlag, vcov = vcov))
  if (class(result)[1] == "try-error") {
    result <- list(rho = matrix(rep(NA, 4),
                                ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                   coefficients = matrix(rep(NA, 18),
                                         ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                   lik = NA,
                   iter = NA)
    if (vcov) {
      result$info = NA
      result$vcov = NA
    }
  }
  return(result)
}


#' @title Inverse digamma function
#' 
#' @description inverse of digamma. digamma function is the first derivative of gamma function divided by gamma function.
#' 
#' @param y a numeric vector.
#' @examples 
#' idigamma(2)
#' plot(digamma, 0.1, 3)
#' plot(idigamma, -10.4, 0.9)
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' 
#' @export
idigamma <- function(y) {
  x = 1
  Vectorize(inv_digamma)(x = x, y = y)
}

# opt <- function(b1, expt, a) {
#   opt_b1(b1, expt, a)
#   return(c(a, b1))
# }


#' @title The bivariate zero-inflated negative binomial distribution 
#'    - Standard error estimation
#' 
#' @description Standard errors of the BZINB distribution parameter estimates are 
#'    calculated based on maximum likelihood estimation. If \code{param} is \code{NULL},
#'    the parameters are first estimated by \code{bzinb} function.
#' 
#' @param xvec,yvec a pair of bzinb random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(a0, a1, a2, b1, b2, p1, p2, p3, p4)}). 
#'    See \code{\link{bzinb}} for more detail.
#' @param a0,a1,a2 shape parameters of the latent gamma variables. They must be positive.
#' @param b1,b2 scale parameters for the latent gamma variables. They must be positive.
#' @param p1,p2,p3,p4 proportions summing up to 1 (\code{p1 + p2 + p3 + p4 = 1}). 
#' \code{p1} is the probability of both latent Poisson variables being observed. 
#' \code{p2} is the probability of only the first Poisson variables being observed.
#' \code{p3} is the probability of only the second Poisson variables being observed, and
#' \code{p4} is the probability of both Poisson variables being dropped out. 
#' @param ... Other arguments passed on to \code{bzinb} function, when \code{param} is 
#'    \code{NULL}.
#' 
#' @return 
#'  Standard error of \code{rho},  \code{logit.rho},  \code{a0, a1, a2, b1, b2, p1, p2, p3},
#'  and \code{p4} estimates, variance-covariance matrix (\code{vcov}) and information matrix.
#'  See \code{\link{bzinb}} for more detail. \code{iter} is \code{NA}, if the \code{param} 
#'  is given.
#'   
#' @examples
#' set.seed(1)
#' data1 <- rbzinb(n = 20, a0 = 1, a1 = 1, a2 = 1, 
#'                 b1 = 1, b2 = 1, p1 = 0.5, p2 = 0.2, 
#'                 p3 = 0.2, p4 = 0.1)
#' bzinb.se(xvec = data1[,1], yvec = data1[,2], 
#'          param = c(5.5, 0.017, 0.017, 0.33, 0.36, 
#'                    0.53, 0.30, 0.08, 0.09))
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references 
#'  Cho, H., Liu, C., Preisser, J., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#' @export
bzinb.se <- function(xvec, yvec, a0, a1, a2, b1, b2, p1, p2, p3, p4, 
                     param = NULL, ...) {
  .check.input(xvec, yvec)
  if (!is.null(param)) {
    if (length(param) != 9) {stop("length(initial) must be 9.")}
    .check.ab(param[1:5])
    .check.p(param[6:9])
  }
  
  if (any(!is.finite(param))) {return(rep(NA, 10))}
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  n <- sum(xy.reduced$freq)
  n.reduced <- as.integer(length(xy.reduced$freq))
  
  if (is.null(param)) {
    warning("param was not provided. It will be estimated.")
    est <- bzinb(xvec, yvec, ...)
    param <- unlist(est$coefficients[,1])
    iter <- est$iter
  } else {
    iter <- NA
  }
  
  if (any(is.na(param))) {
    warning("params include NA.")
    result <- list(rho = matrix(rep(NA, 4),
                                ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                   coefficients = matrix(rep(NA, 18),
                                         ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                   lik = NA, iter = NA, info = NA, vcov = NA)
    return(result)
  }
  
  rho <- param[1]/sqrt((param[1] + param[2]) * (param[1] + param[3])) * 
    sqrt(param[4] *param[5] /(param[4] + 1) /(param[5] + 1))
  logit.rho <- qlogis(rho)
  
  expt = setNames(as.double(rep(0, 12)), expt.names)
  s_i = setNames(as.double(rep(0, 8)), abp.names[-9])
  info <- matrix(0, ncol = 8, nrow = 8, dimnames = list(abp.names[-9], abp.names[-9]))
  dBvZINB_Expt_vec (xvec = xy.reduced$x, yvec = xy.reduced$y, freq = xy.reduced$freq, n = n.reduced, 
                    a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5], 
                    p1 = param[6], p2 = param[7], p3 = param[8], p4 = param[9], 
                    expt = expt, s_i = s_i, info = info, se = 1, bnb = 1)
  
  
  # inverse of info
  qr.info <- try(qr(info))
  if (class(qr.info)[1] == "try-error") {
    warning ("The information matrix has NA/NaN/Inf and thus the standard error is not properly estimatd.")
    std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
    cov.mat <- NA
  } else if (qr(info)$rank < 8) {
    warning ("The information matrix is (essentially) not full rank, and thus the standard error is not reliable.")
    std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
    cov.mat <- NA
  } else {
    cov.mat <- try(solve(info))
    if (class(cov.mat)[1] == "try-error") {
      std.param = setNames(rep(NA, 11), c(abp.names, "rho", "logit.rho"))
      cov.mat <- NA
    } else {
      # variance of p4 hat
      var.p4 <- sum (cov.mat[6:8, 6:8]) # = sum_i,j cov(pi, pj)
      
      # variance of rho hat
      d.g <- rho * c(1/param[1] - 1/{2*(param[1] + param[2])} - 1/{2*(param[1] + param[3])}, 
                     - 1/{2*(param[1] + param[2])}, 
                     - 1/{2*(param[1] + param[3])},
                     1/{2 *param[4] *(param[4] + 1)}, 
                     1/{2 *param[5] *(param[5] + 1)})
      
      var.rho <- t(d.g) %*% cov.mat[1:5, 1:5] %*% d.g
      
      # variance of logit(rho hat)
      var.logit.rho <- var.rho / rho^2 / (1-rho)^2
      # std.param = sqrt(c(setNames(diag(cov.mat), abp.names[1:8]), 
      #                    p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
      std.param = sqrt(c(diag(cov.mat), 
                         p4 = var.p4, rho=var.rho, logit.rho = var.logit.rho))
    }
  } 
  
  result <- list(rho = matrix(c(rho, logit.rho, std.param[c("rho", "logit.rho")]),
                              ncol = 2, dimnames = list(c("rho", "logit.rho"), c("Estimate", "Std.err"))),
                 coefficients = matrix(c(param, std.param[1:9]),
                                       ncol = 2, dimnames = list(abp.names, c("Estimate", "Std.err"))), 
                 lik = expt[1],
                 iter = iter,
                 info = info,
                 vcov = cov.mat)
  return(result)
}
