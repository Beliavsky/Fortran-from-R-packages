\name{stdEr}
\alias{stdEr}
\alias{stdEr.default}
\alias{stdEr.lm}

\title{Standard deviations}
\description{
   Extract standard deviations from estimated models.
}
\usage{
  stdEr(x, ...)
\method{stdEr}{default}(x, \dots)
\method{stdEr}{lm}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, such as created by \code{\link{lm}}}
  \item{\dots}{further arguments for methods}
}
\details{
  \code{stdEr} is a generic function with methods for objects of "lm"
  class.  The default method returns the square root of the diagonal of
  the variance-covariance matrix.
}
\value{
  numeric, the estimated standard errors of the coefficients.  
}
\author{
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link{vcov}}, \code{\link{summary}}.}

\examples{
data(cars)
lmRes <- lm(dist ~ speed, data=cars)
stdEr( lmRes )
}

\keyword{methods}
