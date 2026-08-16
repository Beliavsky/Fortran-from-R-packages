\name{jacobi.g.weight}
\alias{jacobi.g.weight}
\title{ Weight function for the Jacobi polynomial }
\description{
  This function returns the value of the weight function for the order \eqn{k}
  Jacobi polynomial, \eqn{G_k \left( {p,q,x} \right)}.
}
\usage{
jacobi.g.weight(x,p,q)
}
\arguments{
  \item{x}{ the function argument which can be a vector }
  \item{p}{ the first polynomial parameter }
  \item{q}{ the second polynomial parameter }
}
\details{
  The function takes on non-zero values in the interval \eqn{ \left( 0,1 \right) }. The formula
  used to compute the weight function is as follows.
  
  \eqn{w\left( x \right) = \left( {1 - x} \right)^{p - q} \;x^{q - 1}}
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
### compute the Jacobi G weight function for argument values
### between 0 and 1
### parameter p is 3 and q is 2
###
x <- seq( 0, 1, .01 )
y <- jacobi.g.weight( x, 3, 2 )
}
\keyword{ math }
