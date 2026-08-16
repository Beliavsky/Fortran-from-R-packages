\name{poly.calc}
\alias{poly.calc}
\alias{poly.from.zeros}
\alias{poly.from.roots}
\alias{poly.from.values}
\title{
Calculate Polynomials from Zeros or Values
}
\description{
  Calculate either the monic polynomial with specified zeros, or the
  Lagrange interpolation polynomial through the (x,y) points.
}
\usage{
poly.calc(x, y, tol=sqrt(.Machine$double.eps), lab=dimnames(y)[[2]])
}
\arguments{
  \item{x}{numeric vector specifying either the zeros of the desired
    polynomial if this is the only non-missing argument, or the	x-values
    for Lagrange interpolation.}
  \item{y}{numeric vector or matrix specifying the y-values for the
    Lagrange interpolation polynomial.  If \code{y} is a matrix,
    \code{nrow(y)} must equal \code{length(x)}, and each column of
    \code{y} is used separately with \code{x}.}
  \item{tol}{An absolute value tolerance, below which coefficients are
    treated as zero.}
  \item{lab}{If \code{y} is a matrix, lab is used as the names vector
    for the list result.}
}
\details{
  If \code{y} is a matrix, the result is a list of polynomials using
  each column separately.

  If \code{x} only is given, repeated zeros are allowed.  If \code{x}
  and \code{y} are given, repeated values in the \code{x} vector must
  have identical \code{y} values associated with them (up to
  \code{tol}), otherwise the first y-value only is used and a warning
  is issued.
}
\value{
  Either a polynomial object, or a list of polynomials, as appropriate.
  In the latter case the object is of class \code{"polylist"}.
}
\seealso{
  \code{\link{polynomial}}
}
\examples{
poly.calc(rep(1,3))
## -1 + 3*x - 3*x^2 + x^3 
poly.calc(0:4, (0:4)^2 + 1)
## 1 + x^2 
poly.calc(0:4, cbind(0:4, (0:4)^2 + 1), lab = letters[1:2])
## List of polynomials:
## $a:
## x
##
## $b:
## 1 + x^2 
}
\keyword{symbolmath}
