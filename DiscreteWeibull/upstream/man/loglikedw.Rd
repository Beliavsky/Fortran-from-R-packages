\name{loglikedw}
\alias{loglikedw}
\title{
Loglikelihood function
}
\description{
Loglikelihood function (changed in sign) for the type 1 discrete Weibull distribution
}
\usage{
loglikedw(par, x, zero = FALSE)
}
\arguments{
  \item{par}{
the vector of parameters, \eqn{q} and \eqn{\beta}
}
  \item{x}{
the vector of sample values
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
}
\value{
the value of the loglikelihood function (changed in sign) for the observed sample \code{x} under the parameters \code{par} 
}
\author{
Alessandro Barbiero
}
\seealso{\code{\link{estdweibull}}
}
\examples{
x <- c(1,1,1,2,2,2,2,2,2,3,4,4,5,6,8)
-loglikedw(c(0.8, 1), x) # loglikelihood function for q=0.8 and beta=1
-loglikedw(c(0.4, 2), x) # loglikelihood function for q=0.4 and beta=2
par <- estdweibull(x, "ML")# parameter estimates derived by the ML method
par
-loglikedw(par, x) # the maximum value of the loglikelihood function
}
\keyword{distribution}
\keyword{htest}

