##########################################################################################
## 3. Bivariate Negative Binomial (BNB)
##########################################################################################

#' @import Rcpp

### 1. Density, likelihood, gradient :  This function affects bzinb3, bzinb4
dbnb <- function(x, y, a0, a1, a2, b1, b2, log=FALSE) {
  p1 = (b1 + b2 + 1) /(b1 + 1); p2 = (b1 + b2 + 1) /(b2 + 1)
  adj1 <- adj3 <- 0
  l1 <- function(k, m, adj1) exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p1) - adj1)
  l2 <- function(m) sum(sapply(0:x, l1, m = m, adj1 = adj1))
  l2.y <- l2(y)
  while (l2.y > 1e+200) {adj1 = adj1 + 200; l2.y <- l2(y)}
  l3 <- function(m, adj3) log(l2(m)) + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p2) - adj3
  l4 <- function() sum(sapply(0:y, function(s) exp(l3(s, adj3 = adj3))))
  l4.m <- l4()
  while (l4.m > 200) {adj3 = adj3 + 200; l4.m <- l4()}
  # cat("adj1 = ", adj1, "\n"); 
  # cat("adj3 = ", adj3, "\n"); 
  l4.m <- log(l4.m) - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1*log(1+b1) - a2*log(1+b2)  + adj1 + adj3
  return(ifelse(log, l4.m, exp(l4.m)))
}
dbnb.vec <- Vectorize(dbnb)


#' @rdname bnb
#' @name bnb
#' @aliases rbnb
#' @aliases bnb
#' @aliases lik.bnb
#' 
#' @title The bivariate negative binomial distribution
#' 
#' @description random generation (\code{rbnb}), maximum likelihood estimation (\code{bnb}), 
#'    and log-likelihood. (\code{lik.bnb})  for the bivariate negative binomial 
#'    distribution with parameters equal to \code{(a0, a1, a2, b1, b2)}.
#' 
#' @param xvec,yvec a pair of bnb random vectors. nonnegative integer vectors. 
#'    If not integers, they will be rounded to the nearest integers.
#' @param param a vector of parameters (\code{(a0, a1, a2, b1, b2)}). 
#'    Either \code{param} or individual parameters (\code{a0, a1, a2, b1, b2}) 
#'    need to be provided.
#' @param a0,a1,a2 shape parameters of the latent gamma variables. must be positive.
#' @param b1,b2 scale parameters for the latent gamma variables. must be positive.
#' @param n number of observations.
#' @param em if \code{TRUE} in \code{bnb}, EM algorithm is applied. Otherwise, direct opitmation is used.
#' @param tol,maxiter,vcov,initial,showFlag optional arguments applied only when \code{em} is \code{TRUE} in \code{bnb}.
#' 
#' @return 
#'  \itemize{
#'    \item \code{rbnb} gives a pair of random vectors following BNB distribution.
#'    \item \code{bnb} gives the maximum likelihood estimates of a BNB pair. 
#'    Standard error and covariance matrix are provided when \code{em} is \code{TRUE}.
#'    \item \code{lik.bnb} gives the log-likelihood of a set of parameters for a BNB pair.
#'  }
#' 
#' @examples
#' # generating a pair of random vectors
#' set.seed(1)
#' data1 <- rbnb(n = 100, a0 = 2, a1 = 1, a2 = 1, 
#'                 b1 = 1, b2 = 1)
#' 
#' lik.bnb(xvec = data1[, 1], yvec = data1[ ,2], 
#'           a0 = 1, a1 = 1, a2 = 1, b1 = 1, b2 = 1) 
#' 
#' bnb(xvec = data1[,1], yvec = data1[,2], showFlag = FALSE)
#' 
#' @author Hunyong Cho, Chuwen Liu, Jinyoung Park, and Di Wu
#' @references
#'  Cho, H., Liu, C., Preisser, J., and Wu, D. (In preparation), "A bivariate 
#'  zero-inflated negative binomial model for identifying underlying dependence"
#' 
#' @import Rcpp
#' @export
#' @useDynLib bzinb
#'
lik.bnb <- function(xvec, yvec, a0, a1, a2, b1, b2, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 5) stop("length(param) must be 5.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
  }
  sum(log(dbnb.vec(xvec, yvec, a0, a1, a2, b1, b2)))
}



#' @export
#' @rdname bnb
rbnb <- function(n, a0, a1, a2, b1, b2, param = NULL) {
  if (!is.null(param)) {
    if (length(param) != 5) stop("length(param) must be 5.")
    a0 = param[1]; a1 = param[2]; a2 = param[3]; b1 = param[4]; b2 = param[5]
  }
  if (length(n) != 1) {stop("length(n) must be 1.")}
  .check.ab(c(a0, a1, a2, b1, b2))
  
  rmat <- matrix(rgamma(n*3, shape = c(a0, a1, a2), rate = 1/b1), n, 3, byrow=TRUE)
  rmat2 <- rmat
  rmat2[,3] <- rmat[,1] + rmat[,3]
  rmat2[,2] <- rmat[,1] + rmat[,2]
  rmat2 <- rmat2[,2:3]
  rmat2[,2] <- rmat2[,2]*b2/b1
  xy <- matrix(rpois(n*2, rmat2), n, 2)
  colnames(xy) <- c("x", "y")
  
  return(xy)
}

dbnb.gr <- function(x, y, a0, a1, a2, b1, b2) {
  p1 = (b1 + b2 + 1) /(b1 + 1); p2 = (b1 + b2 + 1) /(b2 + 1)
  gr.b1.1 <- x/b1 - (x + y + a0) /(b1 + b2 + 1) - a1 / (b1 + 1)
  gr.b1 <- function(k, m) {(k + m) / (b1 + b2 + 1) - k / (b1 + 1) + gr.b1.1}
  gr.b2.1 <- y/b2 - (x + y + a0) /(b1 + b2 + 1) - a2 / (b2 + 1)
  gr.b2 <- function(k, m) {(k + m) / (b1 + b2 + 1) - m / (b2 + 1) + gr.b2.1}
  gr.a0.1 <- - digamma(a0) - log(1 + b1 + b2)
  gr.a0 <- function(k, m) {if (x==0 & y==0) {- log(1 + b1 + b2)} else  {digamma(x +y +a0 -k -m) + gr.a0.1}}

  gr.a1.1 <- - digamma(a1) - log(1 + b1)
  gr.a1 <- function(k) {if (k==0) {- log(1 + b1)} else  {digamma(a1 + k) + gr.a1.1}}
  gr.a2.1 <- - digamma(a2) - log(1 + b2)
  gr.a2 <- function(m) {if (m==0) {- log(1 + b2)} else  {digamma(a2 + m) + gr.a2.1}}
  
  adj.A <- adj.B1 <- 0
  
  l1 <- function(k, m, adjj = 0) {
    exp(lgamma(a1 + k) - lgamma(k+1) - lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) + k *log(p1)
        + lgamma(m +a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + m *log(p2) - adjj)
  }
  l2 <- - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2)
  if (l2 < - 200) {
    adj.B1 = ((-l2 - 200) %/% 100) * 100 # prevent exp(l1.B) from being 0
    l2 = l2 + adj.B1
  }
  l2 <- exp(l2)
  
  l1.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))
  while (log(sum(l1.mat)) > 250) {
    adj.A = adj.A + 200
    l1.mat <- sapply(0:x, function(k) sapply(0:y, l1, k = k, adjj = adj.A))  # %>% print
  }
  
# cat("line 127\n")
# cat("adjA = ", adj.A, "adjB1 = ", adj.B1, "\n")
# print(l1.mat)
  l1.mat <- (l1.mat * l2)
# cat("line 131\n")
# print(l1.mat)
  l1.sum <- sum(l1.mat) 
  
  gr.b1.mat <- sapply(0:x, function(k) sapply(0:y, gr.b1, k=k))
  gr.b1.mat <- l1.mat * gr.b1.mat
  gr.b1.sum <- sum(gr.b1.mat)/sum(l1.mat)   # (adj.A - adj.B1) cancels out. No need for correction.
  
  gr.b2.mat <- sapply(0:x, function(k) sapply(0:y, gr.b2, k=k))
  gr.b2.mat <- l1.mat * gr.b2.mat
  gr.b2.sum <- sum(gr.b2.mat)/sum(l1.mat)   # (adj.A - adj.B1) cancels out. No need for correction.
    
  gr.a0.mat <- sapply(0:x, function(k) sapply(0:y, gr.a0, k=k))
  gr.a0.mat <- l1.mat * gr.a0.mat
  gr.a0.sum <- sum(gr.a0.mat)/sum(l1.mat)   # (adj.A - adj.B1) cancels out. No need for correction.
  
  gr.a1.mat <- matrix(sapply(0:x, gr.a1),x+1, y+1)
  gr.a1.mat <- l1.mat * t(gr.a1.mat)
  gr.a1.sum <- sum(gr.a1.mat)/sum(l1.mat)   # (adj.A - adj.B1) cancels out. No need for correction.
  
  gr.a2.mat <- matrix(sapply(0:y, gr.a2), y+1, x+1)
  gr.a2.mat <- l1.mat * gr.a2.mat
  gr.a2.sum <- sum(gr.a2.mat)/sum(l1.mat)   # (adj.A - adj.B1) cancels out. No need for correction.
  result <- list(logdensity=log(l1.sum) + adj.A - adj.B1, gradient = c(gr.a0.sum, gr.a1.sum, gr.a2.sum, gr.b1.sum, gr.b2.sum))
  return(result)
}
dbnb.gr.vec <- Vectorize(dbnb.gr)


### 2. MLE
#' @export
#' @rdname bnb
bnb <- function (xvec, yvec, em = TRUE, tol = 1e-8, maxiter = 50000, vcov = TRUE, initial = NULL, showFlag = FALSE) {
  .check.input(xvec, yvec)
  if (!is.null(initial)) {
    if (length(initial) != 5) {stop("length(initial) must be 5.")}
    .check.ab(initial[1:5])
  }
  xvec = as.integer(round(xvec, digits = 0))
  yvec = as.integer(round(yvec, digits = 0))
  
  # initial guess
  xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
  s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
  cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
  initial <- rep(NA,5)
  initial[4] <- s2.x /xbar
  initial[5] <- s2.y /ybar
  initial[2:3] <- c(xbar,ybar)/initial[4:5]
  initial[1] <- min(initial[2:3]) * abs(cor.xy)
  initial[2:3] <-  initial[2:3] - initial[1]
  initial <- pmax(initial, 1e-5)
  
  if (em) {
    initial = c(initial, 1, 0, 0, 0)
    result <- try(bzinb.base(xvec, yvec, initial = initial, tol = tol, maxiter = maxiter, 
                             showFlag = showFlag, vcov = vcov, bnb = 1))
    if (class(result)[1] == "try-error") {
      result <- list(coefficients = matrix(rep(NA, 10),
                                           ncol = 2, dimnames = list(abp.names[1:5], c("Estimate", "Std.err"))), 
                     lik = NA,
                     iter = NA)
      if (vcov) {
        result$info = NA
        result$vcov = NA
      }
    } else {
      result$rho <- NULL
      result$coefficients <- result$coefficients[1:5, ]
    }
  } else {
    result <- bnb.direct (xvec, yvec, tol = tol, showFlag = showFlag)
  }
  return(result)
}
bnb.direct <- function (xvec, yvec, tol = 1e-8, showFlag=FALSE) {
  .check.input(xvec, yvec)
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  if (max(xvec)==0 & max(yvec)==0) return(c(a0 = NA, a1 = NA, a2 = NA, b1 = NA, b2 = NA))
  # tmp.xy.reduced <<- xy.reduced

  fn.1 = function (param) {
    val <- dbnb.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5])
    lik <- sapply(val[1,],cbind) %*% xy.reduced$freq
    return(lik)
  }
  gr.1 = function (param) {
    val <- dbnb.gr.vec( x = xy.reduced$x, y = xy.reduced$y,  a0 = param[1], a1 = param[2], a2 = param[3], b1 = param[4], b2 = param[5])
    lik <- sapply(val[1,],cbind) %*% xy.reduced$freq
    gr <- sapply(val[2,],cbind) %*% xy.reduced$freq
    return(gr)
  }

  
  #log-scaled params: param.l
  fn.log = function (param.l) { fn.1 (exp(param.l))}
  gr.log = function (param.l) { 
    if (showFlag) {message(paste(c("a0", "a1", "a2", "b1", "b2"), round(exp(param.l), 5), collapse = ", "))}
    (as.vector(gr.1 (exp(param.l))) * exp(param.l)) 
  }
  
  # initial guess
  xbar <- mean(xvec); ybar <- mean(yvec); xybar <- mean(c(xbar, ybar))
  s2.x <- var(xvec); s2.y <- var(yvec); if(is.na(s2.x)) {s2.x <- s2.y <- 1}
  cor.xy <- cor(xvec,yvec); if (is.na(cor.xy)) {cor.xy <- 0}
  initial <- rep(NA,5)
  initial[4] <- s2.x /xbar
  initial[5] <- s2.y /ybar
  initial[2:3] <- c(xbar,ybar)/initial[4:5]
  initial[1] <- min(initial[2:3]) * abs(cor.xy)
  initial[2:3] <-  initial[2:3] - initial[1]
  initial <- pmax(initial, 1e-5)
  # tmp.initial <<- initial
  result <- try(exp(optim(par = log(initial), fn = fn.log, gr = gr.log, control=list(fnscale=-1, abstol = tol), method = "BFGS")$par), silent=TRUE)
  if (class(result)[1] == "try-error") {
    initial = rep(1,5)
    result <- exp(optim(par = log(initial), fn = fn.log, gr = gr.log, control=list(fnscale=-1, abstol = tol), method = "BFGS")$par)
  }
  result <- c(a0 = result[1], a1 = result[2], a2 = result[3], b1 = result[4], b2 = result[5])
  return(result)
}


### 3. Deviance
dev.bnb <- function(xvec, yvec, param = NULL, a0 = NULL, a1 = NULL, a2= NULL, b1 = NULL, b2 = NULL) {
  .check.input(xvec, yvec)
  
  # If params = NULL, apply bnb. if else, apply those params
  if (is.null (param)) { 
    if (is.null (a0) | is.null (a1) | is.null (a2) | is.null (b1) | is.null (b2)) {
      param = bnb (xvec = xvec, yvec = yvec)
    }
    else { param = c(a0, a1, a2, b1, b2)}
  }
  # Log-likelihood of the bnb model
  lik.model <- lik.bnb (xvec = xvec, yvec = yvec, param = param)
  
  # Reduced calculation
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$freq <- as.numeric(as.character(xy.reduced$freq))
  
  # Saturated model bzip params
  bnb.vec <- Vectorize(bnb)
  param.sat <- t(bnb.vec(xy.reduced$x, xy.reduced$y))
  param.sat <- do.call(rbind, param.sat)  #"matrix of lists" into "a matrix"
  lik.sat   <- sum(dbnb.vec(x= xy.reduced$x, y = xy.reduced$y, a0 = param.sat[,1], a1 = param.sat[,2], a2 = param.sat[,3], b1 = param.sat[,4], b2 = param.sat[,5],  log = TRUE) * xy.reduced$freq)
  
  return(data.frame(model.likelihood = lik.model, satrtd.likelihood = lik.sat, deviance = 2*(lik.sat - lik.model)))
}