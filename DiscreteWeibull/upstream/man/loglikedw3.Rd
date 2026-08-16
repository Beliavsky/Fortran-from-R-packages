\name{loglikedw3}
\alias{loglikedw3}
\title{
Loglikelihood function
}
\description{
Loglikelihood function (changed in sign) for the type 3 discrete Weibull distribution
}
\usage{
loglikedw3(par, x)
}
\arguments{
  \item{par}{
the vector of parameters, \eqn{c} and \eqn{\beta}
}
  \item{x}{
the vector of sample values
}
}
\value{
the value of the loglikelihood function (changed in sign) for the observed sample \code{x} under the parameters \code{par} 
}
\author{
Alessandro Barbiero
}
\seealso{\code{\link{estdweibull3}}
}
\examples{
n <- 20
c <- 1/3
beta <- 2/3
x <- rdweibull3(n, c, beta)
par <- estdweibull3(x, "ML")
par
-loglikedw3(par, x)
}
\keyword{distribution}
\keyword{htest}

