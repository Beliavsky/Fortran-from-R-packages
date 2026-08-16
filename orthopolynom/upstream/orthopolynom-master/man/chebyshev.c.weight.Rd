\name{chebyshev.c.weight}
\alias{chebyshev.c.weight}
\title{ Weight function for the Chebyshev polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  Chebyshev polynomial of the first kind, \eqn{C_k \left( x \right) }.
}
\usage{
chebyshev.c.weight(x)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -2,2 \right) }.  The formula
  used to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = \frac{1}
{{\sqrt {1 - \frac{{x^2 }}
{4}} }}}
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
### compute the C Chebyshev weight function for arguments between -3 and 3
###
x <- seq( -3, 3, .01 )
y <- chebyshev.c.weight( x )
plot( x, y )
}
\keyword{ math }
