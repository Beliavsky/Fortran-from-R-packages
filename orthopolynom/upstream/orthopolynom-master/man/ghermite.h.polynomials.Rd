\name{ghermite.h.polynomials}
\alias{ghermite.h.polynomials}
\title{ Create list of generalized Hermite polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} generalized Hermite polynomials,   \eqn{H_k^{\left( \mu  \right)} \left( x \right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
ghermite.h.polynomials(n, mu, normalized = FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{mu}{ numeric value for the polynomial parameter }
  \item{normalized}{ boolean value which, if TRUE, returns recurrence relations for normalized polynomials }
}
\details{
  The parameter \eqn{\mu} must be greater than -0.5.
  The function \code{ghermite.h.recurrences} produces a data frame with the recurrence relation parameters
  for the polynomials.  If the \code{normalized} argument is FALSE, the
  function \code{orthogonal.polynomials} is used to construct the list of orthogonal polynomial objects.
  Otherwise, the function \code{orthonormal.polynomials} is used to construct the
  list of orthonormal polynomial objects.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{order 0 generalized Hermite polynomial}
  \item{2 }{order 1 generalized Hermite polynomial}
  ...
  \item{n+1 }{order \eqn{n} generalized Hermite polynomial}
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
\code{\link{ghermite.h.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized generalized Hermite polynomials of orders 0 to 10
### polynomial parameter is 1.0
###
normalized.p.list <- ghermite.h.polynomials( 10, 1, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of unnormalized generalized Hermite polynomials of orders 0 to 10
### polynomial parameter is 1.0
###
unnormalized.p.list <- ghermite.h.polynomials( 10, 1, normalized=FALSE )
print( unnormalized.p.list )
}
\keyword{ math }
