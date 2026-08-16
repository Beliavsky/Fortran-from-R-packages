\name{jacobi.p.inner.products}
\alias{jacobi.p.inner.products}
\title{ Inner products of Jacobi polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} Jacobi polynomial, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
jacobi.p.inner.products(n,alpha,beta)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ numeric value for the first polynomial parameter }
  \item{beta}{ numeric value for the first polynomial parameter }
}
\details{
  The formula used to compute the innser products is as follows.
  
  \eqn{h_n  = \left\langle {P_n^{\left( {\alpha ,\beta } \right)} |P_n^{\left( {\alpha ,\beta } \right)} } \right\rangle  = \frac{{2^{\alpha  + \beta  + 1} }}
{{2\,n + \alpha  + \beta  + 1}}\frac{{\Gamma \left( {n + \alpha  + 1} \right)\,\Gamma \left( {n + \beta  + 1} \right)}}
{{n!\;\Gamma \left( {n + \alpha  + \beta  + 1} \right)}}}.
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
### generate the inner product vector for the P Jacobi polynomials of orders 0 to 10
###
h <- jacobi.p.inner.products( 10, 2, 2 )
print( h )
}
\keyword{ math }
