\name{bvls}
\alias{bvls}
\title{The Stark-Parker implementation of bounded-variable least squares}
\description{
   An R interface to the Stark-Parker implementation of bounded-variable
   least squares that solves the least squares problem
    \eqn{\min{\parallel A x - b \parallel_2}} 
  with the constraint \eqn{l \le x \le u}, where
  \eqn{l,x,u \in R^n, b \in R^m}  and \eqn{A} is an
  \eqn{m \times n} matrix.
}
\usage{
bvls(A, b, bl, bu, key=0, istate=rep(0,ncol(A)+1))
}
\arguments{
 \item{A}{numeric matrix with \code{m} rows and \code{n} columns}
 \item{b}{numeric vector of length \code{m} }
 \item{bl}{numeric vector of length \code{n} specifying the lower bound
   on each element of \code{x}}
 \item{bu}{numeric vector of length \code{n} specifying the upper bound
   on each element of \code{x}}
 \item{key}{If \code{key > 0} the routine initializes using the
    user's guess about
   which components of  \code{x}  are active, i.e. are strictly within their
   bounds, which are at their lower bounds, and which are at their 
   upper bounds.  This information is supplied through the array  
   \code{istate}.} 
 \item{istate}{vector of length \code{ncol(A)+1}.  If \code{key > 0},
   \code{istate} is as follows: the last
   contains the total number of
   components at their bounds (the bound variables).  The absolute
   values of the first \code{nbound <- tail(istate,1)} entries of
   \code{istate} are the indices of these bound components of
   \code{x}.  The sign of \code{istate[1:nbound]} indicates whether
   \code{x(abs(istate[1:nbound]))} is at its upper or lower bound.
   \code{istate[1:nbound]} is positive if the component is at its upper
   bound, negative if the component is at its lower bound.
   \code{istate[(nbound+1):ncol(A)]} contain the indices of the
   components of \code{x} that are active (i.e. are expected to lie
   strictly within their bounds).  When \code{key > 0}, the routine
   initially sets the active components to the averages of their upper
   and lower bounds.
 }
}
\value{
  \code{bvls} returns an object of class \code{"bvls"}.
  
  The generic assessor functions \code{coefficients},
  \code{fitted.values}, \code{deviance} and \code{residuals} extract
  various useful features of the value returned by \code{bvls}.

  An object of class \code{"bvls"} is a list containing the
  following components:

  \item{x}{the parameter estimates.}
  \item{deviance}{the residual sum-of-squares.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted}{the fitted values.}
}
\references{
 Stark PB, Parker RL (1995). Bounded-variable least-squares:
 an algorithm and applications, Computational Statistics, 10, 129-141.
}
\source{
  This is an R interface to the Fortran77 code accompanying the article
  referenced below by Stark PB, Parker RL (1995), and distributed via
  the \bold{statlib} on-line software repository at Carnegie Mellon
  University (URL \url{http://lib.stat.cmu.edu/general/bvls}).  The code
  was modified slightly to allow compatibility with the gfortran
  compiler.  The authors have agreed to distribution under GPL version
  2 or newer.
  
}
\seealso{
  the method \code{"L-BFGS-B"} for \link{optim},
  \link[quadprog]{solve.QP}, \link[nnls]{nnls}
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
X <- matrix(nrow = 50, ncol = 3) 
wavenum <- seq(18000,28000, length=nrow(X))
location <- c(25000, 22000) 
delta <- c(1000,1000)
Xcolfunc <- function(wavenum, location, delta)
  exp( - log(2) * (2 * (wavenum - location)/delta)^2)
for(i in 1:2) X[,i] <- Xcolfunc(wavenum, location[i], delta[i])

X[1:40,3] <- Xcolfunc(wavenum, 23000, 1000)[11:nrow(X)]
X[41:nrow(X),3]<- - Xcolfunc(wavenum, 23000, 1000)[21:30]

## set seed for reproducibility
set.seed(3300)

## simulated data is the product of A and X with some
## spherical Gaussian noise added 
matdat <- A \%*\% t(X) + .005 * rnorm(nrow(A) * nrow(X))

## estimate the rows of X using BVLS criteria 
bvls_sol <- function(matdat, A) {
  X <- matrix(0, nrow = ncol(matdat), ncol = ncol(A) )
  bu <- c(Inf,Inf,.75)
  bl <- c(0,0,-.75)
  for(i in 1:ncol(matdat)) 
     X[i,] <- coef(bvls(A,matdat[,i], bl, bu))
  X
}
X_bvls <- bvls_sol(matdat,A) 

matplot(X,type="p",pch=20)
matplot(X_bvls,type="l",pch=20,add=TRUE)
legend(10, -.5,
c("bound <= zero", "bound <= zero", "bound <= -.75 <= .75"),
col = c(1,2,3), lty=c(1,2,3),
text.col = "blue")

\dontrun{ 
## can solve the same problem with L-BFGS-B algorithm
## but need starting values for x 
bfgs_sol <- function(matdat, A) {
  startval <- rep(0, ncol(A))
  fn1 <- function(par1, b, A) sum( ( b - A \%*\% par1)^2)
  X <- matrix(0, nrow = ncol(matdat), ncol = 3)
  bu <- c(1000,1000,.75)
  bl <- c(0,0,-.75)
  for(i in 1:ncol(matdat))  
    X[i,] <-  optim(startval, fn = fn1, b=matdat[,i], A=A,
                  upper = bu, lower = bl,
                  method="L-BFGS-B")$par
    X
}
X_bfgs <- bfgs_sol(matdat,A) 

## the RMS deviation under BVLS is less than under L-BFGS-B 
sqrt(sum((X - X_bvls)^2)) < sqrt(sum((X - X_bfgs)^2))

## and L-BFGS-B is much slower 
system.time(bvls_sol(matdat,A))
system.time(bfgs_sol(matdat,A))
}
}
\keyword{optimize}
