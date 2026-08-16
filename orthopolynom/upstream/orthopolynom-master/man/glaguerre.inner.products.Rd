\name{glaguerre.inner.products}
\alias{glaguerre.inner.products}
\title{ Inner products of generalized Laguerre polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} generalized Laguerre polynomial, \eqn{L_n^{\left( \alpha \right)} \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
glaguerre.inner.products(n,alpha)
}
\arguments{
  \item{n}{ integer highest polynomial order }
  \item{alpha}{ polynomial parameter }
}
\details{
  The formula used to compute the inner products is as follows.
  
  \eqn{h_n  = \left\langle {L_n^{\left( \alpha  \right)} |L_n^{\left( \alpha  \right)} } \right\rangle  = \frac{{\Gamma \left( {\alpha  + n + 1} \right)}}
{{n!}}}.
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
### generalized Laguerre polynomial inner products of orders 0 to 10
### polynomial parameter is 1.
###
h <- glaguerre.inner.products( 10, 1 )
print( h )
}
\keyword{ math }
