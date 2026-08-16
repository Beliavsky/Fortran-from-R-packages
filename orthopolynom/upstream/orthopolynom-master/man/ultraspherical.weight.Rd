\name{ultraspherical.weight}
\alias{ultraspherical.weight}
\title{ Weight function for the ultraspherical polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} ultraspherical polynomial, \eqn{C_k^{\left( \alpha  \right)} \left( x \right)}.
}
\usage{
ultraspherical.weight(x,alpha)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{alpha}{ polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -1,1 \right) }. The formula used
  to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = \left( {1 - x^2 } \right)^{\alpha  - 0.5} }
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
### compute the ultraspherical weight function for arguments between -2 and 2
### polynomial parameter is 1.0
###
x <- seq( -2, 2, .01 )
y <- ultraspherical.weight( x, 1 )
plot( x, y )
}
\keyword{ math }
