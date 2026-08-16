\name{ghermite.h.recurrences}
\alias{ghermite.h.recurrences}
\title{ Recurrence relations for generalized Hermite polynomials }
\description{
  This function returns a data frame with \eqn{n + 1} rows and four named columns containing
  the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of 
  the recurrence relations for the order \eqn{k} generalized Hermite polynomial, \eqn{H_k^{\left( \mu  \right)} \left( x \right)}, 
  and for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
ghermite.h.recurrences(n, mu, normalized = FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{mu}{ numeric value for the polynomial parameter }
  \item{normalized}{ \code{normalized} boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\details{
  The parameter \eqn{\mu} must be greater than -0.5.
}
\value{
  A data frame with the recurrence relation parameters.
}
\references{
Alvarez-Nordase, R., M. K. Atakishiyeva and N. M. Atakishiyeva, 2004. A q-extension of
the generalized Hermite polynomials with continuous orthogonality property on R,
\emph{International Journal of Pure and Applied Mathematics}, 10(3), 335-347.

Abramowitz, M. and I. A. Stegun, 1968. \emph{Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables}, Dover Publications, Inc., New York.

Courant, R., and D. Hilbert, 1989. \emph{Methods of Mathematical Physics},
John Wiley, New York, NY.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{ghermite.h.inner.products}}
}
\examples{
###
### generate the recurrences data frame for 
### the normalized generalized Hermite polynomials
### of orders 0 to 10.
### polynomial parameter value is 1.0
###
normalized.r <- ghermite.h.recurrences( 10, 1, normalized=TRUE )
print( normalized.r )
###
### generate the recurrences data frame for 
### the unnormalized generalized Hermite polynomials
### of orders 0 to 10.
### polynomial parameter value is 1.0
###
unnormalized.r <- ghermite.h.recurrences( 10, 1, normalized=FALSE )
print( unnormalized.r )
}
\keyword{ math }
