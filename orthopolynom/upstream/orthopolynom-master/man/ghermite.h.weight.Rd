\name{ghermite.h.weight}
\alias{ghermite.h.weight}
\title{ Weight function for the generalized Hermite polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  generalized Hermite polynomial, \eqn{H_k^{\left( \mu  \right)} \left( x \right) }.
}
\usage{
ghermite.h.weight(x, mu)
}
\arguments{
  \item{x}{ a numeric vector function argument }
  \item{mu}{ polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -\infty,\infty \right) }.
  The parameter \eqn{\mu} must be greater than -0.5.  The formula used to compute the
  generalized Hermite weight function is as follows.
  
  \eqn{w\left( {x,\mu } \right) = \left| x \right|^{2\;\mu } \;\exp \left( { - x^2 } \right)}
}
\value{
  The value of the weight function
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
\author{ Frederick Novomestky \email{ fnovomes@poly.edu } }
\examples{
###
### compute the generalized Hermite weight function for argument values 
### between -3 and 3
###
x <- seq( -3, 3, .01 )
y <- ghermite.h.weight( x, 1 )
}
\keyword{ math }
