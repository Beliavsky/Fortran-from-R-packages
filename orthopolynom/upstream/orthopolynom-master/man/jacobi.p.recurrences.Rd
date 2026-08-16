\name{jacobi.p.recurrences}
\alias{jacobi.p.recurrences}
\title{ Recurrence relations for Jacobi polynomials }
\description{
  This function returns a data frame with \eqn{n + 1} rows and four named columns containing 
  the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of 
  the recurrence relations for the order \eqn{k} Jacobi polynomial, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)},
  and for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
jacobi.p.recurrences(n, alpha, beta, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ numeric value for the first polynomial parameter }
  \item{beta}{ numeric value for the second polynomial parameter }
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
\code{\link{jacobi.p.inner.products}},
\code{\link{pochhammer}}
}
\examples{
###
### generate the recurrences data frame for 
### the normalized Jacobi P polynomials
### of orders 0 to 10.
### parameter a is 2 and parameter b is 2
###
normalized.r <- jacobi.p.recurrences( 10, 2, 2, normalized=TRUE )
print( normalized.r )
###
### generate the recurrences data frame for 
### the unnormalized Jacobi P polynomials
### of orders 0 to 10.
### parameter a is 2 and parameter b is 2
###
unnormalized.r <- jacobi.p.recurrences( 10, 2, 2, normalized=FALSE )
print( unnormalized.r )
}
\keyword{ math }
