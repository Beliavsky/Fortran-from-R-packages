\name{orthogonal.polynomials}
\alias{orthogonal.polynomials}
\title{ Create orthogonal polynomials }
\description{
  Create list of orthogonal polynomials from  the following recurrence relations for
  \eqn{k = 0,\;1,\; \ldots ,\;n}.
  
  \deqn{c_k p_{k+1}\left( x \right) = \left( d_k + e_k x \right) p_k \left( x \right) - f_k p_{k-1} \left( x \right)}
  
  We require that \eqn{p_{-1} \left( x \right) = 0} and \eqn{p_0 \left( x \right) = 1}.  
  The coefficients are the column vectors \eqn{{\bf{c}}}, \eqn{{\bf{d}}}, \eqn{{\bf{e}}} and \eqn{{\bf{f}}}.
}
\usage{
orthogonal.polynomials(recurrences)
}
\arguments{
  \item{recurrences}{ a data frame containing the parameters of the orthogonal polynomial recurrence relations }
}
\details{
  The argument is a data frame with \eqn{n + 1} rows and four named columns.  
  The column names are \code{c}, \code{d}, \code{e} and \code{f}.  
  These columns correspond to the column vectors described above.
}
\value{
  A list of \eqn{n + 1} polynomial objects
  \item{1 }{Order 0 orthogonal polynomial}
  \item{2 }{Order 1 orthogonal polynomial}
  ...
  \item{n+1}{Order \eqn{n} orthogonal polynomial}
}
\references{
Abramowitz, M. and I. A. Stegun, 1968. \emph{Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables}, Dover Publications, Inc., New York.

Courant, R., and D. Hilbert, 1989. \emph{Methods of Mathematical Physics},
John Wiley, New York, NY.

Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, 1992. 
\emph{Numerical Recipes in C}, Cambridge University Press, Cambridge, U.K.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate the recurrence relations for T Chebyshev polynomials of orders 0 to 10
###
r <- chebyshev.t.recurrences( 10, normalized=FALSE )
print( r )
###
### generate the list of orthogonal polynomials
###
p.list <- orthogonal.polynomials( r )
print( p.list )
}
\keyword{ math }
