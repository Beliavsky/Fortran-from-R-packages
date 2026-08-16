\name{chebyshev.c.recurrences}
\alias{chebyshev.c.recurrences}
\title{ Recurrence relations for Chebyshev polynomials }
\description{
  This function returns a data frame with \eqn{n + 1} rows and four named columns containing 
  the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of 
  the recurrence relations for the order \eqn{k} Chebyshev polynomial of the first kind, \eqn{C_k \left( x \right)}, 
  and for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
chebyshev.c.recurrences(n, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
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
\code{\link{chebyshev.c.inner.products}}
}
\examples{
###
### generate the recurrences data frame for 
### the normalized Chebyshev C polynomials
### of orders 0 to 10.
###
normalized.r <- chebyshev.c.recurrences( 10, normalized=TRUE )
print( normalized.r )
###
### generate the recurrences data frame for 
### the normalized Chebyshev C polynomials
### of orders 0 to 10.
###
unnormalized.r <- chebyshev.c.recurrences( 10, normalized=FALSE )
print( unnormalized.r )
}
\keyword{ math }
