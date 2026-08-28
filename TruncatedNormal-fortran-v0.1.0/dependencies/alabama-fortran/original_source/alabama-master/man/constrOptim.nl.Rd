\name{constrOptim.nl}
\alias{constrOptim.nl}
\alias{adpbar}
\alias{augpen}
\alias{alabama}
\title{
Nonlinear optimization with constraints
}
\description{
Augmented Lagrangian Adaptive Barrier Minimization Algorithm for optimizing smooth nonlinear objective functions with constraints. Linear or nonlinear equality and inequality constraints are allowed.
}
\usage{
constrOptim.nl(par, fn, gr = NULL, 
hin = NULL, hin.jac = NULL, heq = NULL, heq.jac = NULL, 
control.outer=list(), control.optim = list(), ...)
}
\arguments{
\item{par}{starting vector of parameter values; initial vector must be "feasible"}
\item{fn}{
Nonlinear objective function that is to be optimized. 
      A scalar function that takes a real vector as argument and 
      returns a scalar that is the value of the function at that point 
      (see details).}
  \item{gr}{
The gradient of the objective function \code{fn} evaluated at the 
      argument.  This is a vector-function that takes a real 
      vector as argument and returns a real vector of the same length.  
      It defaults to "NULL", which means that gradient is evaluated numerically.  Computations are dramatically faster in high-dimensional problems when the exact gradient is provided.  See *Example*. 
}
\item{hin}{
a vector function specifying inequality constraints such that hin[j] > 0 for all j
}
\item{hin.jac}{
Jacobian of \code{hin}.  If unspecified, it will be computed using finite-difference, but computations will be faster if specified.
}
\item{heq}{
a vector function specifying equality constraints such that heq[j] = 0 for all j
}
\item{heq.jac}{
Jacobian of \code{heq}.  If unspecified, it will be computed using finite-difference, but computations will be faster if specified.
}
\item{control.outer}{A list of control parameters to be used by the outer loop in \code{constrOptim.nl}.  See *Details* for more information.}
\item{control.optim}{A list of control parameters to be used by the unconstrained optimization algorithm in the inner loop. Identical to that used in \code{optim}.}
\item{...}{Additional arguments passed to \code{fn}, \code{gr}, \code{hin}, \code{heq}.  All of them must accept any specified arguments, either explicitly or by having a \dots argument, but they do not need to use them all.} 
}
\details{
Argument \code{control.outer} is a list specifing any changes to default values of algorithm control parameters for the outer loop.  Note that the names of these must be specified completely.  Partial matching will not work. The list items are as follows:

\code{mu0}: A scaling parameter for barrier penalty for inequality constraints.

\code{sig0}: A scaling parameter for augmented lagrangian for equality constraints

\code{eps}: Tolerance for convergence of outer iterations of the barrier and/or augmented lagrangian algorithm

\code{itmax}: Maximum number of outer iterations.

\code{trace}: A logical variable indicating whether information on outer iterations should be printed out.  If TRUE, at each outer iteration information is displayed on: (i) how well the inequality and equalities are satisfied, (ii) current parameter values, and (iii) current objective function value.

\code{method}: Unconstrained optimization algorithm in optim() to be used; default is the "BFGS" variable metric method.

\code{NMinit}: A logical variable indicating whether "Nelder-Mead" algorithm should be used in optim() for the first outer iteration.
}

\value{
A list with the following components:
\item{par}{Parameters that optimize the nonlinear objective function, satisfying constraints, if convergence is successful.}

\item{value}{The value of the objective function at termination.}

\item{convergence}{An integer code indicating type of convergence.  \code{0} indicates successful convergence. Positive integer codes indicate failure to converge.    
}

\item{message}{Text message indicating the type of convergence or failure.    
}

\item{outer.iterations}{Number of outer iterations} 

\item{lambda}{Value of augmented Lagrangian penalty parameter}

\item{sigma}{Value of augmented Lagrangian penalty parameter for the quadratic term}

\item{barrier.value}{Reduction in the value of the function from its initial value. This is negative in maximization.}

\item{K}{Residual norm of equality constraints.  Must be small at convergence.}

\item{counts}{A vector of length 2 denoting the number of times the objective \code{fn} and the \code{gr} were evaluated, respectively.}
}
\references{

Lange K, \emph{Optimization}, 2004, Springer.

Madsen K, Nielsen HB, Tingleff O, \emph{Optimization With Constraints}, 2004, IMM, Technical University of Denmark.

}
\author{
Ravi Varadhan, Center on Aging and Health, Johns Hopkins University.
}
\seealso{
See Also \code{\link{auglag}}, \code{\link{constrOptim}}.
}

\examples{

fn <- function(x) (x[1] + 3*x[2] + x[3])^2 + 4 * (x[1] - x[2])^2

gr <- function(x) {
g <- rep(NA, 3)
g[1] <- 2*(x[1] + 3*x[2] + x[3]) + 8*(x[1] - x[2]) 
g[2] <- 6*(x[1] + 3*x[2] + x[3]) - 8*(x[1] - x[2]) 
g[3] <- 2*(x[1] + 3*x[2] + x[3])
g
}

heq <- function(x) {
h <- rep(NA, 1)
h[1] <- x[1] + x[2] + x[3] - 1
h
}


heq.jac <- function(x) {
j <- matrix(NA, 1, length(x))
j[1, ] <- c(1, 1, 1)
j
}

hin <- function(x) {
h <- rep(NA, 1)
h[1] <- 6*x[2] + 4*x[3] - x[1]^3 - 3
h[2] <- x[1]
h[3] <- x[2]
h[4] <- x[3]
h
}


hin.jac <- function(x) {
j <- matrix(NA, 4, length(x))
j[1, ] <- c(-3*x[1]^2, 6, 4)
j[2, ] <- c(1, 0, 0)
j[3, ] <- c(0, 1, 0)
j[4, ] <- c(0, 0, 1)
j
}

set.seed(12)
p0 <- runif(3)
ans <- constrOptim.nl(par=p0, fn=fn, gr=gr, heq=heq, heq.jac=heq.jac, hin=hin, hin.jac=hin.jac) 

# Not specifying the gradient and the Jacobians
set.seed(12)
p0 <- runif(3)
ans2 <- constrOptim.nl(par=p0, fn=fn, heq=heq, hin=hin) 

}

\keyword{optimize}
