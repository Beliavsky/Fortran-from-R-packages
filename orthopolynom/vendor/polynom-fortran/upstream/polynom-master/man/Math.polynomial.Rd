\name{Math.polynomial}
\alias{Math.polynomial}
\title{Math Group Methods for Polynomials}
\description{
  Group method for functions in the Math group.
}
\usage{
\method{Math}{polynomial}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"polynomial"}.}
  \item{\dots}{further arguments to be passed to or from methods, such
    as \code{digits}.}
}
\details{
  Most math group functions are disallowed with polynomial arguments.
  The only exceptions are \code{\link{ceiling}}, \code{\link{floor}},
  \code{\link{round}}, \code{\link{trunc}}, and \code{\link{signif}}
  which may be used to transform the coefficients accordingly.
}
\value{
  A polynomial with transformed coefficients.
}
\seealso{
  \code{\link{Ops.polynomial}},
  \code{\link{Summary.polynomial}}.
}
\examples{
op <- options(digits=18)
p <- poly.from.values(1:4, (2:5)^2)
## 1 + 2.00000000000001*x + x^2
p <- round(p)
## 1 + 2*x + x^2
options(op)
}
\keyword{symbolmath}
