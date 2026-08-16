\name{spherical.inner.products}
\alias{spherical.inner.products}
\title{ Inner products of spherical polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product 
  of an order \eqn{k} spherical polynomial, \eqn{P_k \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
spherical.inner.products(n)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
}
\details{
  The formula used to compute the inner products of the spherical orthogonal polynomials
  is the same as that used for the Legendre orthogonal polynomials.
}
\value{
  A vector with \eqn{n + 1} elements
  \item{1 }{inner product of order 0 orthogonal polynomial}
  \item{2 }{inner product of order 1 orthogonal polynomial}
  ...
  \item{n+1 }{inner product of order \eqn{n} orthogonal polynomial}
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
  \code{\link{legendre.inner.products}}
}
\examples{
###
### generate the inner products vector for the spherical polynomals
### of orders 0 to 10.
###
h <- spherical.inner.products( 10 )
print( h )
}
\keyword{ math }
