\name{jacobi.p.polynomials}
\alias{jacobi.p.polynomials}
\title{ Create list of Jacobi polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} Jacobi polynomials, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
jacobi.p.polynomials(n, alpha, beta, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{alpha}{ numeric value for the first polynomial parameter }
  \item{beta}{ numeric value for the second polynomial parameter }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function \code{jacobi.p.recurrences} produces a data frame with the recurrence relation parameters
  for the polynomials.  If the \code{normalized} argument is FALSE, the
  function \code{orthogonal.polynomials} is used to construct the list of orthogonal polynomial objects.
  Otherwise, the function \code{orthonormal.polynomials} is used to construct the
  list of orthonormal polynomial objects.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{order 0 Jacobi polynomial}
  \item{2 }{order 1 Jacobi polynomial}
  ...
  \item{n+1 }{order \eqn{n} Chebyshev polynomial}
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
\code{\link{jacobi.p.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized Jacobi P polynomials of orders 0 to 10
### first parameter value a is 2 and second parameter value b is 2
###
normalized.p.list <- jacobi.p.polynomials( 10, 2, 2, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of unnormalized Jacobi P polynomials of orders 0 to 10
### first parameter value a is 2 and second parameter value b is 2
###
unnormalized.p.list <- jacobi.p.polynomials( 10, 2, 2, normalized=FALSE )
print( unnormalized.p.list )
}
\keyword{ math }
