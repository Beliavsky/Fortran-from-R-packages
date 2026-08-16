\name{hermite.he.inner.products}
\alias{hermite.he.inner.products}
\title{ Inner products of Hermite polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} Hermite polynomial, \eqn{He_k \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
hermite.he.inner.products(n)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
}
\details{
  The formula used to compute the inner products is as follows.
  
  \eqn{h_n  = \left\langle {He_n |He_n } \right\rangle  = \sqrt {2\,\pi } \;n!}.
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
### scaled Hermite polynomials of orders 0 to 10
###
h <- hermite.he.inner.products( 10 )
print( h )
}
\keyword{ math }
