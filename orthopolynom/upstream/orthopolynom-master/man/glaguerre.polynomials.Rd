\name{glaguerre.polynomials}
\alias{glaguerre.polynomials}
\title{ Create list of generalized Laguerre polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{n} generalized Laguerre polynomials, \eqn{L_n^{\left( \alpha  \right)} \left( x \right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
glaguerre.polynomials(n, alpha, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ numeric value for the polynomial parameter }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function \code{glaguerre.recurrences} produces a data frame with the recurrence relation parameters
  for the polynomials.  If the \code{normalized} argument is FALSE, the
  function \code{orthogonal.polynomials} is used to construct the list of orthogonal polynomial objects.
  Otherwise, the function \code{orthonormal.polynomials} is used to construct the
  list of orthonormal polynomial objects.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{order 0 generalized Laguerre polynomial}
  \item{2 }{order 1 generalized Laguerre polynomial}
  ...
  \item{n+1 }{order \eqn{n} generalized Laguerre polynomial}
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
\code{\link{glaguerre.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized generalized Laguerre polynomials of orders 0 to 10
### polynomial parameter is 1.0
###
normalized.p.list <- glaguerre.polynomials( 10, 1, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of unnormalized generalized Laguerre polynomials of orders 0 to 10
### polynomial parameter is 1.0
###
unnormalized.p.list <- glaguerre.polynomials( 10, 1, normalized=FALSE )
}
\keyword{ math }
