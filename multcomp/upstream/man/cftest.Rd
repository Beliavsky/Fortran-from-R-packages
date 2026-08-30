\name{cftest}
\alias{cftest}
\title{ Testing Estimated Coefficients }
\description{
  A convenience function for univariate testing via z- and t-tests of 
  estimated model coefficients
}
\usage{
cftest(model, parm, test = univariate(), ...)
}
\arguments{
  \item{model}{a fitted model.}
  \item{parm}{a vector of parameters to be tested, either a character vector
              of names or an integer.}
  \item{test}{ a function for computing p values, see \code{\link{summary.glht}}.}  
  \item{\dots}{additional arguments passed to \code{\link{summary.glht}}.}
}
\details{
  The usual z- or t-tests are tested without adjusting for 
  multiplicity.
}
\value{
  An object of class \code{summary.glht}.
}
\seealso{
  \code{\link[lmtest]{coeftest}}
}
\examples{

  lmod <- lm(dist ~ speed, data = cars)
  summary(lmod)
  cftest(lmod)
  
}
\keyword{htest}
