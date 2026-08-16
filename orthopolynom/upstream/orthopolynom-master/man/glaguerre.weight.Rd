\name{glaguerre.weight}
\alias{glaguerre.weight}
\title{ Weight function for the generalized Laguerre polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  generalized Laguerre polynomial, \eqn{L_n^{\left( \alpha  \right)} \left( x \right)}.
}
\usage{
glaguerre.weight(x,alpha)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{alpha}{ polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,\infty \right) }.  The formula
  used to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = e^{ - x} \,x^\alpha }
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
### compute the generalized Laguerre weight function for argument values
### between -3 and 3
### polynomial parameter value is 1.0
###
x <- seq( -3, 3, .01 )
y <- glaguerre.weight( x, 1 )
}
\keyword{ math }
