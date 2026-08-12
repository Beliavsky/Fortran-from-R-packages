\name{change_making}
\alias{change_making}
\title{
  Change Making Problem
}
\description{
  Solves the Change Making problem as an integer linear program.
}
\usage{
  change_making(items, value)
}
\arguments{
  \item{items}{vector of integer numbers greater than zero.}
  \item{value}{integer number}
}
\details{
  The Change Making problem attempts to find a minimal combination of
  items that sum up to a given value. If the items are distinct positive
  integers, the solution is unique.

  If the problem is infeasible, i.e. there is no such combination, the
  returned \code{count} is 0.

  The problem is treated as an Integer Linear Program (ILP) and solved
  with the \code{lp} solver in \code{lpSolve}.
}
\value{
  Returns a list with components \code{count}, the number of items used to
  sum up to the value, and \code{solution}, the number of items used
   per item.
}
\references{
  See the Wikipedia article on the "change making problem".
}
\seealso{
  \code{\link{setcover}} 
}
\examples{
items = c(2, 5, 10, 50, 100)
value = 999
change_making(items, value)

\dontrun{
solutions <- numeric(20)
for (m in 1:20) {
    sol <- change_making(items, m)
    solutions[m] <- sol$count
}
solutions
#>  [1] 0 1 0 2 1 3 2 4 3 1 4 2 5 3 2 4 3 5 4 2
}
}
\keyword{ discrete-optimization }
