\name{assignment}
\alias{assignment}
\title{
 Linear Sum Assignment Problem
}
\description{
  Linear (sum) assignment problem, or LSAP.
}
\usage{
assignment(cmat, dir = "min")
}
\arguments{
  \item{cmat}{quadratic (numeric) matrix, the cost matrix.}
  \item{dir}{direction, can be "min" or "max".}
}
\details{
  Solves the linear (sum) assignment problem for quadratic matrices.
  Uses the \code{lp.assign} function from the \code{lpSolve} package,
  that is it solves LSAP as a mixed integer linear programming problem.
}
\value{
  List with components \code{perm}, the permutation that defines the
  minimum solution, \code{min}, the minimum value, and \code{err} is
  always \code{0}, i.e. not used at the moment.
}
\references{
  Burkard, R., M. Dell'Amico, and S. Martello (2009). Assignment Problems.
  Society for Industrial and Applied Mathematics (SIAM).

  Martello, S., and P. Toth (1990). Knapsack Problems: Algorithms and
  Computer Implementations. John Wiley & Sons, Ltd.
}
\note{
  Slower than the Hungarian algorithm in package \code{clue}.
}
\seealso{
  \code{clue::solve_LSAP}
}
\examples{
##  Example similar to clue::solve_LSAP
set.seed(8237)
x <- matrix(sample(1:100), nrow = 10)
y <- assignment(x)
# show permutation and check minimum sum
y$perm                          #   7  6 10  5  8  2  1  4  9  3
y$min                           # 173
z <- cbind(1:10, y$perm)
x[z]                            #  16  9 49  6 17 14  1 44 10  7
y$min == sum(x[z])              # TRUE

\dontrun{
##  Example: minimize sum of distances of complex points
n <- 100
x <- rt(n, df=3) + 1i * rt(n, df=3)
y <- runif(n) + 1i * runif(n)
cmat <- round(outer(x, y, FUN = function(x,y) Mod(x - y)), 2)
system.time(T1 <- assignment(cmat))     # elapsed: 0.003
T1$min / 100                            # 145.75

## Hungarian algorithm in package 'clue'
library("clue")
system.time(T2 <- solve_LSAP(cmat))     # elapsed: 0.014
sum(cmat[cbind(1:n, T2)])               # 145.75
}
}
\keyword{ discrete-optimization }
