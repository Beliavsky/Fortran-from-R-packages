\name{f15.3.1}
\alias{f15.3.1}
\alias{hypergeo_integral}
\title{Hypergeometric function using Euler's integral representation}
\description{
Hypergeometric function using Euler's integral representation, evaluated
using numerical contour integrals.
}
\usage{
f15.3.1(A, B, C, z, h = 0)
}
\arguments{
  \item{A,B,C}{Parameters}
  \item{z}{Primary complex argument}
  \item{h}{specification for the path to be taken; see details}
}
\details{
  Argument \code{h} specifies the path to be taken (the path has to
  avoid the point \eqn{1/z}).  If \code{h} is real and of length 1, the
  path taken comprises two straight lines: one from \eqn{0} to
  \eqn{0.5+hi} and one from \eqn{0.5+hi} to \eqn{1} (if \eqn{h=0} the
  integration is performed over a single segment).

  Otherwise, the integration is performed over \code{length(h)+1}
  segments: \eqn{0} to \code{h[1]}, then \code{h[i]} to \code{h[i+1]}
  for \eqn{1\leq i\leq n-1}{1 <= i <= n-1} and finally \code{h[n]} to 1.

  See examples and notes sections below.
}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\author{Robin K. S. Hankin}
\note{

The Mellin-Barnes form is not yet coded up.

}
\seealso{\code{\link{hypergeo}}}
\examples{


# For |z| <1 the path can be direct:
f15.3.1(2,1,2,-1/2) -2/3

# cf identity 07.23.03.0046.01 of Hypergeometric2F1.pdf with b=1




f <- function(h){f15.3.1(1,2,3, z=2, h=h)}

# Winding number [around 1/z] matters:
f(0.5)
f(c(1-1i, 1+1i, -2i))

# Accuracy isn't too bad; compare numerical to analytical result :
f(0.5) - (-1+1i*pi/2)

}
\keyword{math}
