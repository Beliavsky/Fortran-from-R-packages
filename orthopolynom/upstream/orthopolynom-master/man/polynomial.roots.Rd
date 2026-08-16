\name{polynomial.roots}
\alias{polynomial.roots}
\title{ Create a list of polynomial roots }
\description{
  This function returns a list with \eqn{n} elements containing the roots of 
  the order $k$ monic orthogonal polynomials for orders \eqn{k = 0,\;1,\; \ldots ,\;n} 
  using a data frame with the monic polynomial recurrence parameter vectors 
  \eqn{\bf{a}} and \eqn{\bf{b}}
}
\usage{
polynomial.roots(m.r)
}
\arguments{
  \item{m.r}{ monic recurrence data frame with parameters a and b }
}
\details{
  The parameter \code{m.r} is a data frame with $n$+1 rows and two names columns. 
  The columns which are names \code{a} and \code{b} correspond to the above referenced vectors.
  Function \code{jacobi.matrices} is used to create a list of symmetric, 
  tridiagonal Jacobi matrices from these named columns.  The eigenvalues of 
  the \eqn{k \times k} Jacobi matrix are the roots or zeros of the order $k$ 
  monic orthogonal polynomial.
}
\value{
  A list with \eqn{n} elements each of which is a vector of polynomial roots
  \item{1 }{roots of the order 1 monic polynomial}
  \item{2 }{roots of the order 2 monic polynomial}
  ...
  \item{n }{roots of the order \eqn{n} monic polynomial}
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
  \code{\link{monic.polynomial.recurrences}},
  \code{\link{jacobi.matrices}}
}
\examples{
###
### generate the recurrences data frame for
### the normalized Chebyshev polynomials of
### orders 0 to 10
###
r <- chebyshev.t.recurrences( 10, normalized=TRUE )
###
### construct the corresponding monic polynomial
### recurrences
###
m.r <- monic.polynomial.recurrences( r )
###
### obtain the polynomial roots from the monic polynomial
### recurrences
p.roots <- polynomial.roots( m.r )
}
\keyword{ math }
