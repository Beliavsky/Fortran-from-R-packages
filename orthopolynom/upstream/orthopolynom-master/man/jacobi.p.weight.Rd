\name{jacobi.p.weight}
\alias{jacobi.p.weight}
\title{ Weight function for the Jacobi polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k} 
  Jacobi polynomial, \eqn{P_k^{\left( {\alpha ,\beta } \right)} \left( x \right)}.
}
\usage{
jacobi.p.weight(x,alpha,beta)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{alpha}{ the first polynomial parameter }
  \item{beta}{ the second polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( -1,1 \right) }.  The formula
  used to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = \left( {1 - x} \right)^\alpha  \;\left( {1 + x} \right)^\beta }
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
### compute the Jacobi P weight function for argument values
### between -1 and 1
###
x <- seq( -1, 1, .01 )
y <- jacobi.p.weight( x, 2, 2 )
}
\keyword{ math }
