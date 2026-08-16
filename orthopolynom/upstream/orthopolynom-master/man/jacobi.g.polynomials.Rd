\name{jacobi.g.polynomials}
\alias{jacobi.g.polynomials}
\title{ Create list of Jacobi polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} Jacobi polynomials, \eqn{G_k \left( {p,q,x} \right)}, 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
jacobi.g.polynomials(n, p, q, normalized=FALSE)
}
\arguments{
  \item{n}{ integer value for the highest polynomial order }
  \item{p}{ numeic value for the first polynomial parameter }
  \item{q}{ numeric value for the second polynomial parameter }
  \item{normalized}{ a boolean value which, if TRUE, returns a list of normalized orthogonal polynomials }
}
\details{
  The function \code{jacobi.g.recurrences} produces a data frame with the recurrence relation parameters
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
  \item{n+1 }{order \eqn{n} Jacobi polynomial}
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
\code{\link{jacobi.g.recurrences}},
\code{\link{orthogonal.polynomials}},
\code{\link{orthonormal.polynomials}}
}
\examples{
###
### gemerate a list of normalized Jacobi G polynomials of orders 0 to 10
### first parameter value p is 3 and second parameter value q is 2
###
normalized.p.list <- jacobi.g.polynomials( 10, 3, 2, normalized=TRUE )
print( normalized.p.list )
###
### gemerate a list of normalized Jacobi G polynomials of orders 0 to 10
### first parameter value p is 3 and second parameter value q is 2
###
unnormalized.p.list <- jacobi.g.polynomials( 10, 3, 2, normalized=FALSE )
print( unnormalized.p.list )
}
\keyword{ math }
