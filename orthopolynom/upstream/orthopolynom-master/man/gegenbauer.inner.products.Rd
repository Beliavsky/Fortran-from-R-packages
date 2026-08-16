\name{gegenbauer.inner.products}
\alias{gegenbauer.inner.products}
\title{ Inner products of Gegenbauer polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} Gegenbauer polynomial, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
gegenbauer.inner.products(n,alpha)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ numeric value for the polynomial parameter }
}
\details{
  The formula used to compute the inner products is as follows.
  
  \eqn{h_n  = \left\langle {C_n^{\left( \alpha  \right)} |C_n^{\left( \alpha  \right)} } \right\rangle  = \left\{ {\begin{array}{cc}
   {\frac{{\pi \;2^{1 - 2\,\alpha } \,\Gamma \left( {n + 2\,\alpha } \right)}}
{{n!\;\left( {n + \alpha } \right)\,\left[ {\Gamma \left( \alpha  \right)} \right]^2 }}} & {\alpha  \ne 0}  \\
   {\frac{{2\;\pi }}
{{n^2 }}} & {\alpha  = 0}  \\

 \end{array} } \right.}.
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
  \code{\link{ultraspherical.inner.products}}
}
\examples{
###
### generate the inner products vector for the 
### Gegenbauer polynomials of orders 0 to 10
### the polynomial parameter is 1.0
###
h <- gegenbauer.inner.products( 10, 1 )
print( h )
}
\keyword{ math }
