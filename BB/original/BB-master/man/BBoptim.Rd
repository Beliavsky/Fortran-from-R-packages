\name{BBoptim}
\alias{BBoptim}
\title{Large=Scale Nonlinear Optimization - A Wrapper for spg()}
\description{A strategy using different Barzilai-Borwein steplengths to 
  optimize a nonlinear objective function subject to box constraints.}
\usage{
  BBoptim(par, fn, gr=NULL, method=c(2,3,1), lower=-Inf, upper=Inf, 
  	project=NULL, projectArgs=NULL,
	control=list(), quiet=FALSE, ...) 
  }
\arguments{
    \item{par}{A real vector argument to \code{fn}, indicating the initial 
      guess for the root of the nonliinear system of equations \code{fn}.}

    \item{fn}{Nonlinear objective function that is to be optimized.  A scalar 
      function that takes a real vector as argument and returns a scalar 
      that is the value of the function at that point (see details).}

    \item{gr}{The gradient of the objective function \code{fn} evaluated at the 
      argument.  This is a vector-function that takes a real vector as argument 
      and returns a real vector of the same length.  It defaults to 
      \code{NULL}, which means that gradient is evaluated numerically.
      Computations are dramatically faster in high-dimensional problems when 
      the exact gradient is provided.  See *Example*.} 

    \item{method}{A vector of integers specifying which Barzilai-Borwein 
      steplengths should be used in a consecutive manner.  The methods will 
      be used in the order specified.}

    \item{upper}{An upper bound for box constraints. See \code{spg}}
    \item{lower}{An lower bound for box constraints. See \code{spg}}

    \item{project}{The projection function that takes a point in $R^n$ and 
      projects it onto a region that defines the constraints of the problem.
      This is a vector-function that takes a real vector as argument and 
      returns a real vector of the same length. 
      See \code{spg} for more details.} 

    \item{projectArgs}{list of arguments to \code{project}. See \code{spg()} 
       for more details.} 

    \item{control}{A list of parameters governing the algorithm behaviour.  
      This list is the same as that for \code{spg} (excepting 
      the default for \code{trace}).  See \code{details} for 
      important special features of control parameters.}

    \item{quiet}{logical indicating if messages about convergence success or
       failure should be suppressed} 

    \item{...}{arguments passed fn (via the optimization algorithm).}
}

\value{A list with the same elements as returned by \code{spg}.  One additional 
  element returned is \code{cpar} which contains the control parameter settings
  used to obtain successful convergence, or to obtain the best solution in 
  case of failure.}

\details{
This wrapper is especially useful in problems where (\code{spg} is likely 
to experience convergence difficulties.  When \code{spg()} fails, i.e. 
when \code{convergence > 0} is obtained, a user might attempt various strategies
to find a local optimizer. The function \code{BBoptim} tries the following 
sequential strategy:
\enumerate{
\item Try a different BB steplength.  Since the default is \code{method = 2} 
   for \code{dfsane}, BBoptim wrapper tries \code{method = c(2, 3, 1)}.
\item Try a different non-monotonicity parameter \code{M} for each method, 
  i.e. BBoptim wrapper tries \code{M = c(50, 10)} for each BB steplength.
}

The argument \code{control} defaults to a list with values 
 \code{maxit = 1500, M = c(50, 10), ftol=1.e-10, gtol = 1e-05, maxfeval = 10000, 
  maximize = FALSE, trace = FALSE, triter = 10, eps = 1e-07, checkGrad=NULL}.  
  It is recommended that \code{checkGrad} be set to FALSE for high-dimensional 
  problems, after making sure that the gradient is correctly specified. See 
  \code{spg} for additional details about the default.

If \code{control} is specified as an argument, only values which are different
need to be given in the list. See \code{spg} for more details.
}
\seealso{
  \code{\link{BBsolve}},
  \code{\link{spg}},
  \code{\link{multiStart}}
  \code{\link{optim}}
  \code{\link[numDeriv]{grad}}
   }
   
\examples{
# Use a preset seed so test values are reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

rosbkext <- function(x){
# Extended Rosenbrock function
n <- length(x)
j <- 2 * (1:(n/2))
jm1 <- j - 1
sum(100 * (x[j] - x[jm1]^2)^2 + (1 - x[jm1])^2)
}

p0 <- rnorm(50)
spg(par=p0, fn=rosbkext)
BBoptim(par=p0, fn=rosbkext)

# compare the improvement in convergence when bounds are specified
BBoptim(par=p0, fn=rosbkext, lower=0) 

# identical to spg() with defaults
BBoptim(par=p0, fn=rosbkext, method=3, control=list(M=10, trace=TRUE))  
}
\references{ 
   R Varadhan and PD Gilbert (2009),  BB: An R Package for Solving a Large 
   System of Nonlinear Equations and for Optimizing a High-Dimensional 
   Nonlinear Objective Function, \emph{J. Statistical Software}, 32:4, 
   \url{https://www.jstatsoft.org/v32/i04/}
}

\keyword{multivariate}

