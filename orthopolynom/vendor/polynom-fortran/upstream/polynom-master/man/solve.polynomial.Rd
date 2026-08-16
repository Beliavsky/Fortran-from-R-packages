\name{solve.polynomial}
\alias{solve.polynomial}
\title{Zeros of a Polynomial}
\description{
  Find the zeros, if any, of a given polynomial.
}
\usage{
\method{solve}{polynomial}(a, b, \dots)
}
\arguments{
  \item{a}{A polynomial object for which the zeros are required.}
  \item{b}{a numeric value specifying an additional intercept.  If
    given, the zeros of \code{a - b} are found.}
  \item{\dots}{Not used by this method.}
}
\value{
  A numeric vector, generally complex, of zeros.
}
\details{
  This is a method for the generic function \code{\link{solve}}.
    
  The zeros are found as the eigenvalues of the companion matrix,
  sorted according to their real parts.
}
\seealso{
  \code{\link{polyroot}},
  \code{\link{poly.calc}},
  \code{\link{summary.polynomial}}
}
\examples{
p <- polynomial(6:1)
p
## 6 + 5*x + 4*x^2 + 3*x^3 + 2*x^4 + x^5 
pz <- solve(p)
pz
## [1] -1.49180+0.0000i -0.80579-1.2229i -0.80579+1.2229i
## [4]  0.55169-1.2533i  0.55169+1.2533i
## To retrieve the original polynomial from the zeros:
poly.calc(pz)
## Warning: imaginary parts discarded in coercion
## 6 + 5*x + 4*x^2 + 3*x^3 + 2*x^4 + x^5
}
\keyword{symbolmath}
