\name{laguerre.weight}
\alias{laguerre.weight}
\title{ Weight function for the Laguerre polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  Laguerre polynomial, \eqn{L_n \left( x \right)}.
}
\usage{
laguerre.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,\infty \right) }.  The formula
  used to compute the weight function is as follows.
  \eqn{w\left( x \right) = e^{ - x} } 
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
### compute the Laguerre weight function for argument values
### between 0 and 3
x <- seq( -0, 3, .01 )
y <- laguerre.weight( x )
plot( x, y )
}
\keyword{ math }
