\name{jacobi.g.recurrences}
\alias{jacobi.g.recurrences}
\title{ Recurrence relations for Jacobi polynomials }
\description{
  This function returns a data frame with \eqn{n + 1} rows and four named columns containing 
  the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of 
  the recurrence relations for the order \eqn{k} Jacobi polynomial, \eqn{G_k \left( {p,q,x} \right)}, 
  and for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
jacobi.g.recurrences(n, p, q, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{p}{ numeric value for the first polynomial parameter }
  \item{q}{ numeric value for the second polynomial parameter }
  \item{normalized}{ boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\value{
  A data frame with the recurrence relation parameters.
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
\code{\link{jacobi.g.inner.products}},
\code{\link{pochhammer}}
}
\examples{
###
### generate the recurrences data frame for 
### the normalized Jacobi G polynomials
### of orders 0 to 10.
### parameter p is 3 and parameter q is 2
###
normalized.r <- jacobi.g.recurrences( 10, 3, 2, normalized=TRUE )
print( normalized.r )
###
### generate the recurrences data frame for 
### the normalized Jacobi G polynomials
### of orders 0 to 10.
### parameter p is 3 and parameter q is 2
###
unnormalized.r <- jacobi.g.recurrences( 10, 3, 2, normalized=FALSE )
print( unnormalized.r )
}
\keyword{ math }
