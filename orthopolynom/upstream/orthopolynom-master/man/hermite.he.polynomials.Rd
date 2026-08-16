\name{hermite.he.polynomials}
\alias{hermite.he.polynomials}
\title{ Create list of Hermite polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} Hermite polynomials, \eqn{He_k \left( x \right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
hermite.he.polynomials(n, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for thehighest polynomial order }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function \code{hermite.he.recurrences} produces a data frame with the recurrence relation parameters
  for the polynomials.  If the \code{normalized} argument is FALSE, the
  function \code{orthogonal.polynomials} is used to construct the list of orthogonal polynomial objects.
  Otherwise, the function \code{orthonormal.polynomials} is used to construct the
  list of orthonormal polynomial objects.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{order 0 Hermite polynomial}
  \item{2 }{order 1 Hermite polynomial}
  ...
  \item{n+1 }{order \eqn{n} Hermite polynomial}
}
\references{
Abramowitz, M. and I. A. Stegun, 1968. \emph{Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables}, Dover Publications, Inc., New York.

Courant, R., and D. Hilbert, 1989. \emph{Methods of Mathematical Physics},
John Wiley, New York, NY.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{hermite.he.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized Hermite polynomials of orders 0 to 10
###
normalized.p.list <- hermite.he.polynomials( 10, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of unnormalized Hermite polynomials of orders 0 to 10
###
unnormalized.p.list <- hermite.he.polynomials( 10, normalized=FALSE )
print( unnormalized.p.list )
}
\keyword{ math }
