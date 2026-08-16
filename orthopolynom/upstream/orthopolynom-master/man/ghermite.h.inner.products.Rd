\name{ghermite.h.inner.products}
\alias{ghermite.h.inner.products}
\title{ Inner products of generalized Hermite polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} generalized Hermite polynomial,  \eqn{H_k^{\left( \mu  \right)} \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
ghermite.h.inner.products(n, mu)
}
\arguments{
  \item{n}{ \code{n} integer value for the highest polynomial order }
  \item{mu}{ \code{mu} polynomial parameter }
}
\details{
  The parameter \eqn{\mu} must be greater than -0.5.  The formula used to compute the inner
  products is as follows.
  
  \eqn{h_n \left( \mu  \right) = \left\langle {H_m^{\left( \mu  \right)} |H_n^{\left( \mu  \right)} } \right\rangle  = 2^{2\,n} \,\left[ {\frac{n}
{2}} \right]!\;\Gamma \left( {\left[ {\frac{{n + 1}}
{2}} \right] + \mu  + \frac{1}
{2}} \right)}
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
### generalized Hermite polynomials of orders 0 to 10
### polynomial parameter is 1
###
h <- ghermite.h.inner.products( 10, 1 )
print( h )
}
\keyword{ math }
