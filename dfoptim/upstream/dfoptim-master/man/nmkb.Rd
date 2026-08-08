\name{nmk}
\alias{nmk}
\alias{nmkb}
\title{
Nelder-Mead optimziation algorithm for derivative-free optimization 
}
\description{
An implementation of the Nelder-Mead algorithm for derivative-free optimization.  This allows bounds to be placed on parameters. Bounds are enforced by means of a parameter transformation.}
\usage{
nmk(par, fn, control = list(), ...)

nmkb(par, fn, lower=-Inf, upper=Inf, control = list(), ...)
}
\arguments{
\item{par}{A starting vector of parameter values. Must be feasible, i.e. lie strictly between lower and upper bounds.}
\item{fn}{
Nonlinear objective function that is to be optimized. 
      A scalar function that takes a real vector as argument and 
      returns a scalar that is the value of the function at that point 
      (see details).}
\item{lower}{Lower bounds on the parameters.  A vector of the same length as the parameters.  If a single value is specified, it is assumed that the same lower bound applies to all parameters.} 
\item{upper}{Upper bounds on the parameters.  A vector of the same length as the parameters.  If a single value is specified, it is assumed that the same upper bound applies to all parameters.} 
  \item{control}{A list of control parameters.  See *Details* for more information.
}
  \item{\dots}{Additional arguments passed to \code{fn}
}
}
\details{
Argument \code{control} is a list specifing any changes to default values of algorithm control parameters for the outer loop.  Note that the names of these must be specified completely.  Partial matching will not work.  The list items are as follows:

\code{tol} Convergence tolerance.  Iteration is terminated when the absolute difference in function value between successive iteration is below \code{tol}.  Default is 1.e-06.

\code{maxfeval}: Maximum number of objective function evaluations allowed.  Default is min(5000, max(1500, 20*length(par)^2)).

\code{regsimp} A logical variable indicating whether the starting parameter configuration is a regular simplex.  Default is TRUE.

\code{maximize} A logical variable indicating whether the objective function should be maximized.  Default is FALSE.

\code{restarts.max} Maximum number of times the algorithm should be restarted before declaring failure. Default is 3.

\code{trace} A logical variable indicating whether the starting parameter configuration is a regular simplex.  Default is FALSE.
}
\value{
A list with the following components:
\item{par}{Best estimate of the parameter vector found by the algorithm.}

\item{value}{The value of the objective function at termination.}

\item{feval}{The number of times the objective \code{fn} was evaluated.
}

\item{restarts}{The number of times the algorithm had to be restarted when it stagnated.
}

\item{convergence}{An integer code indicating type of convergence.  \code{0} indicates successful convergence. Positive integer codes indicate failure to converge.    
}

\item{message}{Text message indicating the type of convergence or failure.  
}

}
\references{
C.T. Kelley (1999), Iterative Methods for Optimization, SIAM.
}
\author{
Ravi Varadhan <rvaradhan@jhmi.edu>, Johns Hopkins University
URL:http://www.jhsph.edu/agingandhealth/People/Faculty_personal_pages/Varadhan.html
}
\note{
This algorithm is based on the Matlab code of Prof. C.T. Kelley, given in his book "Iterative methods for optimization".
It is implemented here with the permission of Prof. Kelley and SIAM.  However, there are some non-trivial modifications of the algorithm.
}


\seealso{
\code{\link{optim}}, \code{\link{hjk}}, \code{\link{mads}}
}
\examples{
 rosbkext <- function(x){
# Extended Rosenbrock function
 n <- length(x)
 sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
 }

np <- 10
set.seed(123)

p0 <- rnorm(np)
xm1 <- nmk(fn=rosbkext, par=p0) # maximum `fevals' is not sufficient to find correct minimum
xm1b <- nmkb(fn=rosbkext, par=p0, lower=-2, upper=2)

### A non-smooth problem
hald <- function(x) {
#Hald J & Madsen K (1981), Combined LP and quasi-Newton methods 
#for minimax optimization, Mathematical Programming, 20, p.42-62.
	i <- 1:21
	t <- -1 + (i - 1)/10
	f <- (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
	max(abs(f))
	}

p0 <- runif(5)
xm2 <- nmk(fn=hald, par=p0)
xm2b <- nmkb(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)

## Another non-smooth functions
  nsf <- function(x) {
	f1 <- x[1]^2 + x[2]^2
	f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
	f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
	max(f1, f2, f3)
  }
par0 <- c(1, 1)                                 # true min 7.2 at (1.2, 2.4)
nmk(par0, nsf) # fmin=8 at xmin=(2,2)

}
\keyword{optimize}
