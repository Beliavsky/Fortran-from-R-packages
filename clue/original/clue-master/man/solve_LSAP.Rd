\name{solve_LSAP}
\encoding{UTF-8}
\alias{solve_LSAP}
\title{Solve Linear Sum Assignment Problem}
\description{
  Solve the linear sum assignment problem using the Hungarian method.
}
\usage{
solve_LSAP(x, maximum = FALSE)
}
\arguments{
  \item{x}{a matrix with nonnegative entries and at least as many
    columns as rows.}
  \item{maximum}{a logical indicating whether to minimize of maximize
    the sum of assigned costs.}
}
\details{
  If \eqn{nr} and \eqn{nc} are the numbers of rows and columns of
  \code{x}, \code{solve_LSAP} finds an optimal \emph{assignment} of rows
  to columns, i.e., a one-to-one map \code{p} of the numbers from 1 to
  \eqn{nr} to the numbers from 1 to \eqn{nc} (a permutation of these
  numbers in case \code{x} is a square matrix) such that
  \eqn{\sum_{i=1}^{nr} x[i, p[i]]} is minimized or maximized.

  This assignment can be found using a linear program (and package
  \pkg{lpSolve} provides a function \code{lp.assign} for doing so), but
  typically more efficiently and provably in polynomial time
  \eqn{O(n^3)} using primal-dual methods such as the so-called Hungarian
  method (see the references).
}
\value{
  An object of class \code{"solve_LSAP"} with the optimal assignment of
  rows to columns.
}
\references{
  \bibshow{Papadimitriou+Steiglitz:1982}
}
\author{
  Walter Böhm \email{Walter.Boehm@wu.ac.at} kindly provided C code
  implementing the Hungarian method.
}
\examples{
x <- matrix(c(5, 1, 4, 3, 5, 2, 2, 4, 4), nrow = 3)
solve_LSAP(x)
solve_LSAP(x, maximum = TRUE)
## To get the optimal value (for now):
y <- solve_LSAP(x)
sum(x[cbind(seq_along(y), y)])
}
\keyword{optimize}
