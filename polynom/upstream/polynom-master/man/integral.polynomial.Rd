\name{integral.polynomial}
\alias{integral}
\alias{integral.polynomial}
\title{Integrate a Polynomial}
\description{
  Find the integral of a univariate polynomial.
}
\usage{
\method{integral}{polynomial}(expr, limits = NULL, ...)
}
\arguments{
  \item{expr}{an object of class \code{"polynomial"}.}
  \item{limits}{numeric vector of length 2 giving the integration
    limits.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
  If \code{limits} is not given, the integral of \code{p} from 0 to
  \sQuote{x}'.  Otherwise, the integral with the given integration
  limits.
}
\seealso{
  \code{\link{deriv.polynomial}}
}
\examples{
p <- poly.calc(1:5)
p
## -120 + 274*x - 225*x^2 + 85*x^3 - 15*x^4 + x^5
deriv(p)
## 274 - 450*x + 255*x^2 - 60*x^3 + 5*x^4
integral(deriv(p)) - 120
## -120 + 274*x - 225*x^2 + 85*x^3 - 15*x^4 + x^5
}
\keyword{symbolmath}
