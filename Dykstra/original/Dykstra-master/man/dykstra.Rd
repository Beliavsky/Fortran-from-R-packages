\name{dykstra}
\alias{dykstra}

\title{
Solve a Quadratic Programming Problem via Dykstra's Algorithm
}
\description{
This function uses Dykstra's cyclic projection algorithm to solve quadratic programming problems of the form \deqn{- d^T x + (1/2) x^T D x } subject to \eqn{A^T x >= b} where \eqn{D} is a positive definite (or positive semidefinite) matrix. 
}
\usage{
dykstra(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
        maxit = NULL, eps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dmat}{
Quadratic program matrix \eqn{D} of order \eqn{n \times n}.
}
  \item{dvec}{
Quadratic program vector \eqn{d} of length \eqn{n}.
}
  \item{Amat}{
Constraint matrix \eqn{A} of order \eqn{n \times r}.
}
  \item{bvec}{
Constraint vector \eqn{b} of length \eqn{r}. Defaults to vector of zeros.
}
  \item{meq}{
First \code{meq} constraints are equality constraints (remaining are inequality constraints). Defaults to zero.
}
  \item{factorized}{
If \code{TRUE}, argument \code{Dmat} is \eqn{R^{-1}} where \eqn{R^T R = D}.  
}
  \item{maxit}{
Maximum number of iterations (cycles). Defaults to \eqn{30n}.
}
  \item{eps}{
Numeric tolerance. Defaults to \code{n * .Machine$double.eps}.
}
}
\details{
Arguments 1-6 of the \code{\link{dykstra}} function are inspired by (and identical to) the corresponding arguments of the \code{\link[quadprog]{solve.QP}} function in the \bold{quadprog} package.
}
\value{
\item{solution}{Vector \eqn{x} that minimizes quadratic function subject to constraints.}
\item{value}{Value of quadratic function at \code{solution}. Will be \code{NA} if \code{factorized = TRUE}.}
\item{unconstrained}{Vector \eqn{x_0 = D^{-1} d} that minimizes quadratic function ignoring constraints.}
\item{iterations}{Number of iterations (cycles) of the algorithm.}
\item{converged}{\code{TRUE} if algorithm converged. \code{FALSE} if iteration limit exceeded.}
}
\references{
Dykstra, Richard L. (1983). An algorithm for restricted least squares regression. \emph{Journal of the American Statistical Association}, Volume 78, Issue 384, 837-842. doi: 10.1080/01621459.1983.10477029
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
For positive semidefinite \eqn{D}, a small constant is added to each eigenvalue of \eqn{D} before solving the quadratic programming problem.
}
\examples{
###  EXAMPLE 1: Generic Quadratic Programming Problem  ###

# constraint 1 (equality): coefficients sum to 1
# constraints 2-4 (inequality): coefficients non-negative

# define QP problem
Dmat <- diag(3)
dvec <- c(1, 1.5, 1)
Amat <- cbind(rep(1, 3), diag(3))
bvec <- c(1, 0, 0, 0)

# solve QP problem
dykstra(Dmat, dvec, Amat, bvec, meq = 1)

# solve QP problem (factorized = TRUE)
dykstra(Dmat, dvec, Amat, bvec, meq = 1, factorized = TRUE)



###  EXAMPLE 2: Regression with Non-Negative Coefficients  ###

# generate regression data
set.seed(1)
nobs <- 100
nvar <- 5
X <- matrix(rnorm(nobs*nvar), nobs, nvar)
beta <- c(0, 1, 0.3, 0.7, 0.1)
y <- X \%*\% beta + rnorm(nobs)

# define QP problem
Dmat <- crossprod(X)
dvec <- crossprod(X, y)
Amat <- diag(nvar)

# solve QP problem
dykstra(Dmat, dvec, Amat)

# solve QP problem (factorized = TRUE)
Rmat <- chol(Dmat)
Rinv <- solve(Rmat)
dykstra(Rinv, dvec, Amat, factorized = TRUE)



###  EXAMPLE 3: Isotonic Regression  ###

# generate regression data
set.seed(1)
n <- 50
x <- 1:n
y <- log(x) + rnorm(n)

# define QP problem
Dmat <- diag(n)
Amat <- Dmat[, 2:n] - Dmat[, 1:(n-1)]

# solve QP problem
dyk <- dykstra(Dmat, y, Amat)
dyk

# plot results
plot(x, y)
lines(x, dyk$solution)



###  EX 4: Large Non-Negative Quadratic Program  ###

# define QP problem
set.seed(1)
n <- 1000
Dmat <- Amat <- diag(n)
dvec <- runif(n, min = -2)

# solve QP problem with dykstra
dyk <- dykstra(Dmat, dvec, Amat)
dyk

}
\keyword{ optimize }
