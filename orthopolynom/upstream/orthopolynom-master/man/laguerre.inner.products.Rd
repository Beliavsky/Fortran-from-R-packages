\name{laguerre.inner.products}
\alias{laguerre.inner.products}
\title{ Inner products of Laguerre polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} Laguerre polynomial, \eqn{L_n \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
laguerre.inner.products(n)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
}
\details{
  The formula used to compute the inner products is as follows.
  
  \eqn{h_n  = \left\langle {L_n |L_n } \right\rangle  = 1}.
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
\examples{
###
### generate the inner products vector for the
### Laguerre polynomial inner products of orders 0 to 10
###
h <- laguerre.inner.products( 10 )
print( h )
}
\keyword{ math }
