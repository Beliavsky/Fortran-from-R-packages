\name{monic}
\alias{monic}
\title{Monic Polynomials}
\description{
  Convert a polynomial to monic form by dividing by the leading
  coefficient.
}
\usage{
monic(p)
}
\arguments{
  \item{p}{A polynomial.  A warning is issued if the polynomial is
    identically zero.}
}
\details{
  Similar in effect to \code{p/as.numeric(p[length(p)])} but with some
  safeguards against leading zero coefficients.
}
\value{
  A polynomial proportional to \code{p} with leading coefficient 1.
}
\seealso{
  \code{\link{change.origin}}
}
\keyword{symbolmath}
