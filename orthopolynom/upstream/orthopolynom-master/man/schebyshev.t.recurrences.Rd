\name{schebyshev.t.recurrences}
\alias{schebyshev.t.recurrences}
\title{ Recurrence relations for shifted Chebyshev polynomials }
\description{
  This function returns a data frame with \eqn{n + 1} rows and four named columns containing 
  the coefficient vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} of the recurrence relations 
  for the order \eqn{k} shifted Chebyshev polynomial of the first kind, \eqn{T_k^* \left( x \right)}, and 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
schebyshev.t.recurrences(n, normalized)
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

Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, 1992. 
\emph{Numerical Recipes in C}, Cambridge University Press, Cambridge, U.K.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
\code{\link{schebyshev.t.inner.products}}
}
\examples{
###
### generate the recurrence relations for 
### the normalized shifted T Chebyshev polynomials
### of orders 0 to 10
###
normalized.r <- schebyshev.t.recurrences( 10, normalized=TRUE )
print( normalized.r )
###
### generate the recurrence relations for 
### the unnormalized shifted T Chebyshev polynomials
### of orders 0 to 10
###
unnormalized.r <- schebyshev.t.recurrences( 10, normalized=FALSE )
print( unnormalized.r )
}
\keyword{ math }
