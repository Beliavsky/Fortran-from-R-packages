\name{jacobi.g.inner.products}
\alias{jacobi.g.inner.products}
\title{ Inner products of Jacobi polynomials }
\description{
  This function returns a vector with \eqn{n + 1} elements containing the inner product of 
  an order \eqn{k} Jacobi polynomial, \eqn{G_k \left( {p,q,x} \right)}, 
  with itself (i.e. the norm squared) for orders \eqn{k = 0,\;1,\; \ldots ,\;n }.
}
\usage{
jacobi.g.inner.products(n,p,q)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{p}{ numeric value for the first polynomial parameter }
  \item{q}{ numeric value for the first polynomial parameter }
}
\details{
  The formula used to compute the inner products is as follows.
  
  \eqn{h_n  = \left\langle {G_n |G_n } \right\rangle  = \frac{{n!\;\Gamma \left( {n + q} \right)\,\Gamma \left( {n + p} \right)\,\Gamma \left( {n + p - q + 1} \right)}}
{{\left( {2\,n + p} \right)\,\left[ {\Gamma \left( {2\,n + p} \right)} \right]^2 }}}.
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
### G Jacobi polynomials of orders 0 to 10
### parameter p is 3 and parameter q is 2
###
h <- jacobi.g.inner.products( 10, 3, 2 )
print( h )
}
\keyword{ math }
