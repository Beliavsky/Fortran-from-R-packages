\name{hermite.he.weight}
\alias{hermite.he.weight}
\title{ Weight function for the Hermite polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  Hermite polynomial, \eqn{He_k \left( x \right)}.
}
\usage{
hermite.he.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -\infty,\infty \right) }. The
  formula used to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = \exp \left( { - \frac{{x^2 }}{2}} \right)}
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
### compute the scaled Hermite weight function for argument values
### between -3 and 3
###
x <- seq( -3, 3, .01 )
y <- hermite.he.weight( x )
}
\keyword{ math }
