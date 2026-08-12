\name{setcover}
\alias{setcover}
\title{
  Set cover problem
}
\description{
  Solves the Set Cover problem as an integer linear program.
}
\usage{
  setcover(Sets, weights)
}
\arguments{
  \item{Sets}{matrix of 0s and 1s, each line defining a subset.}
  \item{weights}{numerical weights for each subset.}
}
\details{
  The Set Cover problems attempts to find in subsets (of a 'universe')
  a minimal set of subsets that still covers the whole set.

  Each line of the matrix \code{Sets} defines a characteristic function of
  a subset. It is required that each element of the universe is contained
  in at least one of these subsets.

  The problem is treated as an Integer Linear Program (ILP) and solved
  with the \code{lp} solver in \code{lpSolve}.
}
\value{
  Returns a list with components \code{sets}, giving the indices of subsets,
  and \code{objective}, the sum of weights of subsets present in the solution.
}
\references{
  See the Wikipedia article on the "set cover problem".
}
\seealso{
  \code{\link{knapsack}} 
}
\examples{
# Define 12 subsets of universe {1, ..., 10}.
set.seed(7*11*13)
A <- matrix(sample(c(0,1), prob = c(0.8,0.2), size = 120, replace =TRUE),
            nrow = 12, ncol = 10)
sol <- setcover(Sets = A, weights = rep(1, 12))
sol
## $sets
## [1]  1  2  9 12
## $no.sets
##[1] 4

# all universe elements are covered:
colSums(A[sol$sets, ])
## [1] 1 1 2 1 1 1 2 1 1 2
}
\keyword{ discrete-optimization }
