\name{Ops.polynomial}
\alias{Ops.polynomial}
\title{Arithmetic Ops Group Methods for Polynomials}
\description{
  Allows arithmetic operators to be used for polynomial calculations,
  such as addition, multiplication, division, etc.
}
\usage{
\method{Ops}{polynomial}(e1, e2)
}
\arguments{
  \item{e1}{an object of class \code{"polynomial"}.}
  \item{e2}{an object of class \code{"polynomial"}.}
}
\value{
  A polynomial got by performing the operation on the two arguments.
}
\seealso{
  \code{\link{Math.polynomial}},
  \code{\link{Summary.polynomial}}.
}
\examples{
p <- polynomial(c(1, 2, 1))
## 1 + 2*x + x^2
r <- poly.calc(-1 : 1)
## -1*x + x^3
(r - 2 * p)^2
## 4 + 20*x + 33*x^2 + 16*x^3 - 6*x^4 - 4*x^5 + x^6
}
\keyword{symbolmath}
