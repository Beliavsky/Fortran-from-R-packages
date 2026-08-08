\name{multiStart}
\alias{multiStart}
\title{Nonlinear Optimization or Root-Finding with Multiple Starting Values}
\description{Start \code{BBsolve} or \code{BBoptim} from multiple starting 
points to obtain multiple solutions and to test sensitivity to starting values.}
\usage{
  multiStart(par, fn, gr=NULL, action = c("solve", "optimize"), 
	method=c(2,3,1),  lower=-Inf, upper=Inf,
	project=NULL, projectArgs=NULL, 
	control=list(),  quiet=FALSE, details=FALSE, ...) 
  }
\arguments{
    \item{par}{A real matrix, each row of which is an argument to \code{fn},
      indicating initial guesses for solving a nonlinear system \code{fn = 0} 
      or for optimizing the objective function \code{fn}.}

    \item{fn}{see \code{BBsolve} or \code{BBoptim}.} 

    \item{gr}{Only required for optimization.  See \code{BBoptim}.} 

    \item{action}{A character string indicating whether to solve a nonlinear 
      system or to optimize.  Default is ``solve''.} 

    \item{method}{see \code{BBsolve} or \code{BBoptim}.} 

    \item{upper}{An upper bound for box constraints. See \code{spg}}
    \item{lower}{An lower bound for box constraints. See \code{spg}}

    \item{project}{A projection
      function or character string indicating its name. The projection
      function that takes a point in \eqn{R^n}{R^n} and 
      projects it onto a region that defines the constraints of the problem.  
      This is a vector-function that takes a real vector as argument and 
      returns a real vector of the same length. 
      See \code{spg} for more details.} 

    \item{projectArgs}{A list with arguments to the \code{project}  
      function.}

    \item{control}{See \code{BBsolve} and \code{BBoptim}.}

    \item{quiet}{A logical variable (TRUE/FALSE). If \code{TRUE} warnings 
    and some additional information printing are suppressed.  
    Default is \code{quiet = FALSE}
    Note that the \code{control} variable \code{trace} and \code{quiet} 
    affect different printing, so if \code{trace} is not set to \code{FALSE} 
    there will be considerable printed output.}
    \item{details}{Logical indicating if the result should include the full
      result from \code{BBsolve} or \code{BBoptim} for each starting value.}
    \item{...}{arguments passed fn (via the optimization algorithm).}
}
\value{list with elements \code{par}, \code{values}, and \code{converged}. 
  It optionally returns an attribute called ``details'', which is a list as long as 
  the number of starting values, which contains the complete object returned 
  by \code{dfsane} or \code{spg} for each starting value.}

\details{
The optimization or root-finder is run with each row of \code{par} indicating 
initial guesses.
}
\seealso{
  \code{\link{BBsolve}},
  \code{\link{BBoptim}},
  \code{\link{dfsane}},
  \code{\link{spg}}
   }
   
\examples{
# Use a preset seed so the example is reproducable. 
require("setRNG")
old.seed <- setRNG(list(kind="Mersenne-Twister", normal.kind="Inversion",
    seed=1234))

# Finding multiple roots of a nonlinear system
brownlin <- function(x) {
# Brown's almost linear system(A.P. Morgan, ACM 1983)
# two distinct solutions if n is even
# three distinct solutions if n is odd  
  	n <- length(x)
  	f <- rep(NA, n)
	nm1 <- 1:(n-1)
	f[nm1] <- x[nm1] + sum(x) - (n+1)
	f[n] <- prod(x) - 1 
	f
}

p <- 9
n <- 50
p0 <- matrix(rnorm(n*p), n, p)  # n starting values, each of length p
ans <- multiStart(par=p0, fn=brownlin)
pmat <- ans$par[ans$conv, 1:p] # selecting only converged solutions
ord1 <- order(abs(pmat[,1]))
round(pmat[ord1, ], 3)  # all 3 roots can be seen

# An optimization example
rosbkext <- function(x){
n <- length(x)
j <- 2 * (1:(n/2))
jm1 <- j - 1
sum(100 * (x[j] - x[jm1]^2)^2 + (1 - x[jm1])^2)
}

p0 <- rnorm(50)
spg(par=p0, fn=rosbkext)
BBoptim(par=p0, fn=rosbkext)

pmat <- matrix(rnorm(100), 20, 5)  # 20 starting values each of length 5 
ans <- multiStart(par=pmat, fn=rosbkext, action="optimize")
ans
attr(ans, "details")[[1]]  # 

pmat <- ans$par[ans$conv, 1:5] # selecting only converged solutions
round(pmat, 3)
}
\references{ 
   R Varadhan and PD Gilbert (2009),  BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function, \emph{J. Statistical Software}, 32:4, \url{https://www.jstatsoft.org/v32/i04/}
}

\keyword{multivariate}
