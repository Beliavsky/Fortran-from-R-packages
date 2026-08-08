\name{hjk}
\alias{hjk}
\alias{hjkb}
\title{
  Hooke-Jeeves derivative-free minimization algorithm
}
\description{
  An implementation of the Hooke-Jeeves algorithm for derivative-free
  optimization.  A bounded and an unbounded version are provided.
}
\usage{
hjk(par, fn, control = list(), ...)

hjkb(par, fn, lower = -Inf, upper = Inf, control = list(), ...)
}
\arguments{
  \item{par}{Starting vector of parameter values.    The initial vector may lie on the boundary. If \code{lower[i]=upper[i]}
  for some \code{i}, the \code{i}-th component of the solution vector will
  simply be kept fixed.}
  \item{fn}{Nonlinear objective function that is to be optimized. 
            A scalar function that takes a real vector as argument and 
            returns a scalar that is the value of the function at that point.}
  \item{lower, upper}{Lower and upper bounds on the parameters.
                      A vector of the same length as the parameters.
                      If a single value is specified, it is assumed that the 
                      same bound applies to all parameters. The
                      starting parameter values must lie within the bounds.}
  \item{control}{A list of control parameters.
                 See \bold{Details} for more information.}
  \item{\dots}{Additional arguments passed to \code{fn}.}
}
\details{
  Argument \code{control} is a list specifing changes to default values of
  algorithm control parameters.
  Note that parameter names may be abbreviated as long as they are unique.

  The list items are as follows:
  \describe{
  \item{\code{tol}}{Convergence tolerance. Iteration is terminated when the 
  step length of the main loop becomes smaller than \code{tol}. This does
  \emph{not} imply that the optimum is found with the same accuracy.
  Default is 1.e-06.}

  \item{\code{maxfeval}}{Maximum number of objective function evaluations 
  allowed. Default is Inf, that is no restriction at all.}

  \item{\code{maximize}}{A logical indicating whether the objective function 
  is to be maximized (TRUE) or minimized (FALSE). Default is FALSE.}

  \item{\code{target}}{A real number restricting the absolute function value. 
  The procedure stops if this value is exceeded.
  Default is Inf, that is no restriction.}

  \item{\code{info}}{A logical variable indicating whether the step number, 
  number of function calls, best function value, and the first component of
  the solution vector will be printed to the console. Default is FALSE.}
  }
  If the minimization process threatens to go into an infinite loop, set
  either \code{maxfeval} or \code{target}.
}
\value{
  A list with the following components:
  \item{par}{Best estimate of the parameter vector found by the algorithm.}

  \item{value}{value of the objective function at termination.}

  \item{convergence}{indicates convergence (\code{=0}) or not (\code{=1}).}

  \item{feval}{number of times the objective \code{fn} was evaluated.}

  \item{niter}{number of iterations in the main loop.}
}
\references{
  C.T. Kelley (1999), Iterative Methods for Optimization, SIAM.

  Quarteroni, Sacco, and Saleri (2007), Numerical Mathematics, Springer.
}
\author{
  Hans W Borchers  <hwborchers@googlemail.com>
}
\note{
  This algorithm is based on the Matlab code of Prof. C. T. Kelley, given
  in his book ``Iterative methods for optimization".
  It is implemented here with the permission of Prof. Kelley.

  This version does not (yet) implement a cache for storing function values
  that have already been computed as searching the cache makes it slower.
}
\seealso{
  \code{\link{optim}}, \code{\link{nmk}}
}
\examples{
##  Hooke-Jeeves solves high-dim. Rosenbrock function
  rosenbrock <- function(x){
    n <- length(x)
    sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
  }
par0 <- rep(0, 10)
hjk(par0, rosenbrock)

hjkb(c(0, 0, 0), rosenbrock, upper = 0.5)
# $par
# [1] 0.50000000 0.25742722 0.06626892


##  Hooke-Jeeves does not work well on non-smooth functions
  nsf <- function(x) {
	f1 <- x[1]^2 + x[2]^2
	f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
	f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
	max(f1, f2, f3)
  }
par0 <- c(1, 1)                                 # true min 7.2 at (1.2, 2.4)
hjk(par0, nsf) # fmin=8 at xmin=(2,2)
}
\keyword{ optimize }
