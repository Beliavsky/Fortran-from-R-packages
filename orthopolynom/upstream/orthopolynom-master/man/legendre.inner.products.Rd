\name{legendre.inner.products}
\alias{legendre.inner.products}
\title{ Inner products of Legendre polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product 
  of an order \eqn{k} Legendre polynomial, \eqn{P_k \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
legendre.inner.products(n)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
}
\details{
  The formula used compute the inner products is as follows.
  
  \eqn{h_n  = \left\langle {P_n |P_n } \right\rangle  = \frac{2}
{{2\,n + 1}}}. 
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
  \code{\link{spherical.inner.products}}
}
\examples{
###
### compute the inner product for the
###  Legendre polynomials of orders 0 to 1
###
h <- legendre.inner.products( 10 )
print( h )
}
\keyword{ math }
