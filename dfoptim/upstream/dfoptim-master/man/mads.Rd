\name{mads}
\alias{mads}
\alias{mads}
\title{
Mesh Adaptive Direct Searches (MADS) algorithm for derivative-free and black-box optimization 
}
\description{
An implementation of the Mesh Adaptive Direct Searches (MADS) algorithm for derivative-free and black-box optimization. It uses a series of variable size meshes to search the space and to converge to (local) minima with mathematical proof of convergence. It is usable on unbounded and bounded unconstrained problems. The objective function can return ``NA'' if out-of-bound or violating constraints (strict barrier approach for constraints), or a penalty can be added to the objective function.}
\usage{
mads(par, fn, lower=-Inf, upper=Inf, scale=1, control = list(), ...)
}
\arguments{
\item{par}{A starting vector of parameter values. Must be feasible, i.e. lie strictly between lower and upper bounds.}
\item{fn}{
Noisy, non-differentiable, non-convex, piecewise or nonlinear objective function that is to be optimized. It takes a real vector as argument and returns a scalar or ``NA'' that is the value of the function at that point (see details).}
\item{lower}{Lower bounds on the parameters.  A vector of the same length as the parameters.  If a single value is specified, it is assumed that the same lower bound applies to all parameters. If all lower bounds are -Inf and all upper bounds are Inf, then the problem is treated as unbounded.} 
\item{upper}{Upper bounds on the parameters.  A vector of the same length as the parameters.  If a single value is specified, it is assumed that the same upper bound applies to all parameters. If all lower bounds are -Inf and all upper bounds are Inf, then the problem is treated as unbounded.}
\item{scale}{Optional scaling, default is 1.  A vector of the same length as the parameters. If a single value is specified, it is assumed that the same scale factor applies to all parameters.  This scale factor can be customized for each parameter allowing non-proportional moves in the space (normally used for unbounded problems). }  
  \item{control}{A list of control parameters.  See *Details* for more information.  }
  \item{\dots}{Additional arguments passed to \code{fn}
}
}
\details{
Argument \code{control} is a list specifing any changes to default values of algorithm control parameters for the outer loop.  The list items are as follows:

\code{tol} Convergence tolerance.  Iteration is terminated when the absolute difference in function value between successive iteration is below \code{tol}.  Default is 1.e-06.

\code{maxfeval}: Maximum number of objective function evaluations allowed.  Default is 10000).

\code{trace} A logical variable indicating whether information is printed on the console during execution.  Default is TRUE.

\code{maximize} A logical variable indicating whether the objective function should be maximized.  Default is FALSE (hence default is minimization).

\code{pollStyle} A string variable indicating density of the poll set, or, number of vectors in the positive basis. Choices are: ``lite'' (n+1 points) or ``full'' (2n points).  Default is ``lite''.

\code{deltaInit} A numerical value specifying the initial mesh size, between ``tol'' and 1 (mesh size is limited to 1).  Default is 0.01.

\code{expand} A numerical value >1 specifying the expansion (is success) and contraction (if no success) factor of the mesh at the end of an iteration. Default is 4.

\code{lineSearch} A integer value indicating the maximum of search steps to consider. Line search is performed at the end of a successful poll set evaluation, along the line going from last to new ``best'' solution. Stepsize will be automatically increased according to the Fibonacci series. Default is 20. Set to -1 to disable the feature.

\code{seed} Seed value for the internal pseudo random numbers generator. Default is 1138.


}
\value{
A list with the following components:
\item{par}{Best estimate of the parameter vector found by the algorithm.}

\item{value}{The value of the objective function at termination.}

\item{feval}{The number of times the objective \code{fn} was evaluated.
}

\item{convergence}{Final mesh size, should be <tol if successfule convergence. If feval reached maxfeval, then the algorithm did not converge.
}

\item{iterlog}{A dataframe used to log properties of the ``best'' solution at the end of each iteration. 
}

}
\references{
C. Audet and J. E. Dennis, Jr. Mesh adaptive direct search algorithms for constrained optimization. SIAM Journal on Optimization, 17(1): 188-217, 2006.
}
\author{
Vincent Bechard <vincent.bechard@hec.ca>, HEC Montreal (Montreal University)
URL:https://www.linkedin.com/in/vincentbechard 
}
\note{
This algorithm is based on the Lower Triangular method described in the reference.
}


\seealso{
\code{\link{optim}}, \code{\link{hjk}}, \code{\link{nmk}}
}

\examples{
 rosbkext <- function(x){
# Extended Rosenbrock function
 n <- length(x)
 sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
 }

np <- 10
p0 <- rnorm(np)
ans1 <- mads(fn=rosbkext, par=p0, lower=-10, upper=10, scale=1, control=list(trace=FALSE))

### A non-smooth problem from Hock & Schittkowski #78
hs78 <- function(x){
  f <- rep(NA, 3)
  f[1] <- sum(x^2) - 10
  f[2] <- x[2]*x[3] - 5*x[4]*x[5]
  f[3] <- x[1]^3 + x[2]^3 + 1
  F <- prod(x) + 10*sum(abs(f))
  return(F)
}

p0 <- c(-2,1.5,2,-1,-1)
ans2 <- mads(p0, hs78, control=list(trace=FALSE)) #minimum value around -2.81


}
\keyword{optimize}
