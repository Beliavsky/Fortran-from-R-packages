\name{nnnpls}
\alias{nnnpls}
\title{An implementation of least squares with non-negative and non-positive
  constraints}
\description{
  An implementation of an algorithm for linear least squares problems
  with non-negative and non-positive
  constraints based on the Lawson-Hanson
  NNLS algorithm.   Solves \eqn{\min{\parallel A x - b \parallel_2}}
  with the constraint \eqn{x_i \ge 0}
  if \eqn{con_i \ge 0} and \eqn{x_i \le 0} otherwise,  where
  \eqn{x, con \in R^n, b \in R^m}, and \eqn{A} is an \eqn{m \times n} matrix.  
}
\usage{
nnnpls(A, b, con)
}
\arguments{
 \item{A}{numeric matrix with \code{m} rows and \code{n} columns}
 \item{b}{numeric vector of length \code{m} }
 \item{con}{numeric vector of length \code{m} where element \code{i}
   is negative if and only if element \code{i} of the solution vector
   \code{x} should be constrained to non-positive, as opposed to
   non-negative, values. }
} 
\value{
  \code{nnnpls} returns
  an object of class \code{"nnnpls"}.
  
  The generic accessor functions \code{coefficients},
  \code{fitted.values}, \code{deviance} and \code{residuals} extract
  various useful features of the value returned by \code{nnnpls}.

  An object of class \code{"nnnpls"} is a list containing the
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
  This is an R interface to Fortran77 code distributed
  with the book referenced below by Lawson CL, Hanson RJ (1995),
  obtained from Netlib (file \file{lawson-hanson/all}), with some
  trivial modifications to allow for the combination of constraints to
  non-negative and non-positive values, and to return the variable
  NSETP.
}
\seealso{
\link{nnls}, the method \code{"L-BFGS-B"} for \link{optim},
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
## 2 of the Gaussian shapes are non-negative and 1 is non-positive 
X <- matrix(nrow = 51, ncol = 3)
wavenum <- seq(18000,28000, by=200)
location <- c(25000, 22000, 20000) 
delta <- c(3000,3000,3000)
Xcolfunc <- function(wavenum, location, delta)
  exp( - log(2) * (2 * (wavenum - location)/delta)^2)
for(i in 1:3) X[,i] <- Xcolfunc(wavenum, location[i], delta[i])
X[,2] <- -X[,2]

## set seed for reproducibility
set.seed(3300)

## simulated data is the product of A and X with some
## spherical Gaussian noise added 
matdat <- A \%*\% t(X) + .005 * rnorm(nrow(A) * nrow(X))

## estimate the rows of X using NNNPLS criteria 
nnnpls_sol <- function(matdat, A) {
  X <- matrix(0, nrow = 51, ncol = 3)
  for(i in 1:ncol(matdat)) 
     X[i,] <- coef(nnnpls(A,matdat[,i],con=c(1,-1,1)))
  X
}
X_nnnpls <- nnnpls_sol(matdat,A) 

\dontrun{ 
## can solve the same problem with L-BFGS-B algorithm
## but need starting values for x and 
## impose a very low/high bound where none is desired
bfgs_sol <- function(matdat, A) {
  startval <- rep(0, ncol(A))
  fn1 <- function(par1, b, A) sum( ( b - A \%*\% par1)^2)
  X <- matrix(0, nrow = 51, ncol = 3)
  for(i in 1:ncol(matdat))  
    X[i,] <-  optim(startval, fn = fn1, b=matdat[,i], A=A,
              lower=rep(0, -1000, 0), upper=c(1000,0,1000),
              method="L-BFGS-B")$par
    X
}
X_bfgs <- bfgs_sol(matdat,A) 

## the RMS deviation under NNNPLS is less than under L-BFGS-B 
sqrt(sum((X - X_nnnpls)^2)) < sqrt(sum((X - X_bfgs)^2))

## and L-BFGS-B is much slower 
system.time(nnnpls_sol(matdat,A))
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
#  Amat <- diag(c(1,-1,1))
#  for(i in 1:ncol(matdat)) { 
#    dvec <- crossprod(A,matdat[,i]) 
#    X[i,] <- solve.QP(dvec = dvec, bvec = bvec, Dmat=Dmat,
#                      Amat=Amat)$solution
#  }
#  X
# }
# X_quadprog <- quadprog_sol(matdat,A) 

## the RMS deviation under NNNPLS is about the same as under quadprog 
# sqrt(sum((X - X_nnnpls)^2))
# sqrt(sum((X - X_quadprog)^2))

## and quadprog requires about the same amount of time 
# system.time(nnnpls_sol(matdat,A))
# system.time(quadprog_sol(matdat,A))
}
}
\keyword{optimize}
