\name{ultraspherical.polynomials}
\alias{ultraspherical.polynomials}
\title{ Create list of ultraspherical polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} ultraspherical polynomials, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
ultraspherical.polynomials(n, alpha, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ polynomial parameter }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function \code{ultraspherical.recurrences} produces a data frame with the recurrence relation parameters
  for the polynomials.  If the \code{normalized} argument is FALSE, the
  function \code{orthogonal.polynomials} is used to construct the list of orthogonal polynomial objects.
  Otherwise, the function \code{orthonormal.polynomials} is used to construct the
  list of orthonormal polynomial objects.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{order 0 ultraspherical polynomial}
  \item{2 }{order 1 ultraspherical polynomial}
  ...
  \item{n+1 }{order \eqn{n} ultraspherical polynomial}
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
\code{\link{gegenbauer.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized ultra spherical polynomials 
### of orders 0 to 10
###
normalized.p.list <- ultraspherical.polynomials( 10, 1, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of unnormalized ultra spherical polynomials 
### of orders 0 to 10
###
unnormalized.p.list <- ultraspherical.polynomials( 10, 1, normalized=FALSE )
print( unnormalized.p.list )
}
\keyword{ math }
