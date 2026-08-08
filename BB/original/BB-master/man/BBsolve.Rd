\name{BBsolve}
\alias{BBsolve}
\title{Solving Nonlinear System of Equations - A Wrapper for dfsane()}
\description{A strategy using different Barzilai-Borwein steplengths to 
  solve a nonlinear system of equations.}
\usage{
  BBsolve(par, fn, method=c(2,3,1), 
	control=list(), quiet=FALSE, ...) 
  }
\arguments{
    \item{par}{A real vector argument to \code{fn}, indicating the initial guess
       for the root of the nonliinear system of equations \code{fn}.}

    \item{fn}{Nonlinear system of equation that is to be solved. 
      A vector function that takes a real vector as argument and 
      returns a real vector of the same length.}

    \item{method}{A vector of integers specifying which Barzilai-Borwein 
      steplengths should be used in a consecutive manner.  The methods will 
      be used in the order specified.}

    \item{control}{A list of parameters governing the algorithm behaviour.  
      This list is the same as that for \code{dfsane} and \code{sane} (excepting 
      the default for \code{trace}).  
      See \code{details} for important special features of control parameters.}

    \item{quiet}{logical indicating if messages about convergence success or
       failure should be suppressed} 

    \item{...}{arguments passed fn (via the optimization algorithm).}
}

\value{A list with the same elements as returned by \code{dfsane} 
  or \code{sane}.  One additional element returned is \code{cpar} which 
  contains the control parameter settings used to obtain successful 
  convergence, or to obtain the best solution in case of failure. }

\details{
This wrapper is especially useful in problems where the algorithms
(\code{dfsane} or \code{sane}) are likely to experience difficulties in 
convergence.  When these algorithms with default parameters fail, i.e. 
when \code{convergence > 0} is obtained, a user might attempt various 
strategies to find a root of the nonlinear system. The function \code{BBsolve} 
tries the following sequential strategy:
\enumerate{
\item Try a different BB steplength.  Since the default is \code{method = 2} 
  for \code{dfsane}, the BBsolve wrapper tries \code{method = c(2, 1, 3)}. 
\item Try a different non-monotonicity parameter \code{M} for each method, 
  i.e. BBsolve wrapper tries \code{M = c(50, 10)} for each BB steplength.
\item Try with Nelder-Mead initialization.  Since the default for 
  \code{dfsane} is \code{NM = FALSE}, BBsolve does \code{NM = c(TRUE, FALSE)}.
}

The argument \code{control} defaults to a list with values 
 \code{maxit = 1500, M = c(50, 10), tol = 1e-07, trace = FALSE, 
        triter = 10, noimp = 100, NM = c(TRUE, FALSE)}. 
If \code{control} is specified as an argument, only values which are different
need to be given in the list. See \code{dfsane} for more details.
}
\seealso{
  \code{\link{BBoptim}},
  \code{\link{dfsane}},
  \code{\link{sane}}
  \code{\link{multiStart}}
   }
   
\examples{
# Use a preset seed so test values are reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

broydt <- function(x) {
n <- length(x)
f <- rep(NA, n)
h <- 2
f[1] <- ((3 - h*x[1]) * x[1]) - 2*x[2] + 1
tnm1 <- 2:(n-1)
f[tnm1] <- ((3 - h*x[tnm1]) * x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
f[n] <- ((3 - h*x[n]) * x[n]) - x[n-1] + 1
f
}

p0 <- rnorm(50)
BBsolve(par=p0, fn=broydt)  # this works 
dfsane(par=p0, fn=broydt) # but this is highly unliikely to work.

# this implements the 3 BB steplengths with M = 50, and without Nelder-Mead initialization
BBsolve(par=p0, fn=broydt, control=list(M=50, NM=FALSE))

# this implements BB steplength 1 with M = 50 and 10, and both with and 
#   without Nelder-Mead initialization  
BBsolve(par=p0, fn=broydt, method=1, control=list(M=c(50, 10))) 

# identical to dfsane() with defaults
BBsolve(par=p0, fn=broydt, method=2, control=list(M=10, NM=FALSE)) 
}
\references{ 
   R Varadhan and PD Gilbert (2009),  BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function, \emph{J. Statistical Software}, 32:4, \url{https://www.jstatsoft.org/v32/i04/}
}

\keyword{multivariate}
