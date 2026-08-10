\name{nnls}
\alias{nnls}
\title{The Lawson-Hanson NNLS implemention of non-negative least squares}
\description{
  An R interface to the Lawson-Hanson
  NNLS implementation of an algorithm
  for non-negative linear least squares 
  that solves 
  \eqn{\min{\parallel A x - b \parallel_2}} with the
    constraint \eqn{x \ge 0},  where
  \eqn{x \in R^n, b \in R^m}  and \eqn{A} is an \eqn{m \times n} matrix. 
}
\usage{
nnls(A, b)
}
\arguments{
 \item{A}{numeric matrix with \code{m} rows and \code{n} columns}
 \item{b}{numeric vector of length \code{m} }
} 
\value{
 \code{nnls} returns
  an object of class \code{"nnls"}.
  
  The generic accessor functions \code{coefficients},
  \code{fitted.values}, \code{deviance} and \code{residuals} extract
  various useful features of the value returned by \code{nnls}.

  An object of class \code{"nnls"} is a list containing the
  following components:

  
   \item{x}{the parameter estimates.}
  \item{deviance}{the residual sum-of-squares.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted}{the fitted values.}
  \item{mode}{a character vector containing a message regarding why
    termination occured.}
  \item{passive}{vector of the indices of \code{x} that are not bound
    at zero. }
  \item{bound}{vector of the indices of \code{x} that are bound
    at zero.}
  \item{nsetp}{the number of elements of \code{x} that are not bound
  at zero. }
}
\references{
Lawson CL, Hanson RJ (1974). Solving Least Squares Problems. Prentice
Hall, Englewood Cliffs, NJ.

Lawson CL, Hanson RJ (1995). Solving Least Squares Problems. Classics
in Applied Mathematics. SIAM, Philadelphia.
}

\source{
  This is an R interface to the Fortran77 code distributed
  with the book referenced below by Lawson CL, Hanson RJ (1995),
  obtained from Netlib (file \file{lawson-hanson/all}), 
  with a trivial modification to return the variable
  NSETP.
}
\seealso{\link{nnnpls}, the method \code{"L-BFGS-B"} for \link{optim},
  \link[quadprog]{solve.QP}, \link[bvls]{bvls}
} 
\examples{
## simulate a matrix A
## with 3 columns, each containing an exponential decay 
t <- seq(0, 2, by = .04)
k <- c(.5, .6, 1)
A <- matrix(nrow = 51, ncol = 3)
Acolfunc <- function(k, t) exp(-k*t)
for(i in 1:3) A[,i] <- Acolfunc(k[i],t)

## simulate a matrix X
## with 3 columns, each containing a Gaussian shape 
## the Gaussian shapes are non-negative
X <- matrix(nrow = 51, ncol = 3)
wavenum <- seq(18000,28000, by=200)
location <- c(25000, 22000, 20000) 
delta <- c(3000,3000,3000)
Xcolfunc <- function(wavenum, location, delta)
  exp( - log(2) * (2 * (wavenum - location)/delta)^2)
for(i in 1:3) X[,i] <- Xcolfunc(wavenum, location[i], delta[i])

## set seed for reproducibility
set.seed(3300)

## simulated data is the product of A and X with some
## spherical Gaussian noise added 
matdat <- A \%*\% t(X) + .005 * rnorm(nrow(A) * nrow(X))

## estimate the rows of X using NNLS criteria 
nnls_sol <- function(matdat, A) {
  X <- matrix(0, nrow = 51, ncol = 3)
  for(i in 1:ncol(matdat)) 
     X[i,] <- coef(nnls(A,matdat[,i]))
  X
}
X_nnls <- nnls_sol(matdat,A) 

matplot(X_nnls,type="b",pch=20)
abline(0,0, col=grey(.6))

\dontrun{
## can solve the same problem with L-BFGS-B algorithm
## but need starting values for x 
bfgs_sol <- function(matdat, A) {
  startval <- rep(0, ncol(A))
  fn1 <- function(par1, b, A) sum( ( b - A \%*\% par1)^2)
  X <- matrix(0, nrow = 51, ncol = 3)
  for(i in 1:ncol(matdat))  
    X[i,] <-  optim(startval, fn = fn1, b=matdat[,i], A=A,
                   lower = rep(0,3), method="L-BFGS-B")$par
    X
}
X_bfgs <- bfgs_sol(matdat,A) 

## the RMS deviation under NNLS is less than under L-BFGS-B 
sqrt(sum((X - X_nnls)^2)) < sqrt(sum((X - X_bfgs)^2))

## and L-BFGS-B is much slower 
system.time(nnls_sol(matdat,A))
system.time(bfgs_sol(matdat,A))

## can also solve the same problem by reformulating it as a
## quadratic program (this requires the quadprog package; if you
## have quadprog installed, uncomment lines below starting with
## only 1 "#" )

# library(quadprog)

# quadprog_sol <- function(matdat, A) {
#  X <- matrix(0, nrow = 51, ncol = 3)
#  bvec <- rep(0, ncol(A))
#  Dmat <- crossprod(A,A)
#  Amat <- diag(ncol(A))
#  for(i in 1:ncol(matdat)) { 
#    dvec <- crossprod(A,matdat[,i]) 
#    X[i,] <- solve.QP(dvec = dvec, bvec = bvec, Dmat=Dmat,
#                      Amat=Amat)$solution
#  }
#  X
# }
# X_quadprog <- quadprog_sol(matdat,A) 

## the RMS deviation under NNLS is about the same as under quadprog 
# sqrt(sum((X - X_nnls)^2))
# sqrt(sum((X - X_quadprog)^2))

## and quadprog requires about the same amount of time 
# system.time(nnls_sol(matdat,A))
# system.time(quadprog_sol(matdat,A))

}

}
\keyword{optimize}
