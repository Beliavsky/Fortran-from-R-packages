\name{schebyshev.t.polynomials}
\alias{schebyshev.t.polynomials}
\title{ Create list of shifted Chebyshev polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} shifted Chebyshev polynomials of the first kind, \eqn{T_k^* \left( x\right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
schebyshev.t.polynomials(n, normalized)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function \code{schebyshev.t.recurrences} produces a data frame with the recurrence relation parameters
  for the polynomials.  If the \code{normalized} argument is FALSE, the
  function \code{orthogonal.polynomials} is used to construct the list of orthogonal polynomial objects.
  Otherwise, the function \code{orthonormal.polynomials} is used to construct the
  list of orthonormal polynomial objects.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{order 0 shifted Chebyshev polynomial}
  \item{2 }{order 1 shifted Chebyshev polynomial}
  ...
  \item{n+1 }{order \eqn{n} shifted Chebyshev polynomial}
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
\code{\link{schebyshev.u.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized shifted T Chebyshev polynomials of orders 0 to 10
###
normalized.p.list <- schebyshev.t.polynomials( 10, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of unnormalized shifted T Chebyshev polynomials of orders 0 to 10
###
unnormalized.p.list <- schebyshev.t.polynomials( 10, normalized=FALSE )
print( unnormalized.p.list )
}
\keyword{ math }
