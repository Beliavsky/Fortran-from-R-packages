\name{ultraspherical.inner.products}
\alias{ultraspherical.inner.products}
\title{ Inner products of ultraspherical polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} ultraspherical polynomial, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
ultraspherical.inner.products(n,alpha)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ numeric value for the polynomial parameter }
}
\details{
  This function uses the same formula as the function \code{gegenbauer.inner.products}.
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
Formulas, Graphs, and Mathematical Tables}, Dover Publications, Inc., NY.

Courant, R., and D. Hilbert, 1989. \emph{Methods of Mathematical Physics},
John Wiley, New York, NY.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{gegenbauer.inner.products}}
}
\examples{
###
### generate the inner products vector for the
### ultraspherical polynomials of orders 0 to 10.
### the polynomial parameter is 1.0
###
h <- ultraspherical.inner.products( 10, 1 )
print( h )
}
\keyword{ math }
