\name{minuslogRoy}
\alias{minuslogRoy}
\title{
Log-likelihood function
}
\description{
Log-likelihood function (with minus sign) for Roy's bivariate geometric model
}
\usage{
minuslogRoy(x, y, theta1 = 0.5, theta2 = 0.5, theta3 = 1)
}
\arguments{
  \item{x}{
a vector of observed values (non-negative integers)
}
  \item{y}{
a vector of observed values (non-negative integers) of the same length as \code{x}
}
  \item{theta1}{
paramater \eqn{\theta_1}
}
  \item{theta2}{
paramater \eqn{\theta_2}
}
  \item{theta3}{
paramater \eqn{\theta_3}
}
}
\value{
The value of the log-likelihood function, changed in sign
}
\note{
Just to be used inside the \code{estbivgeomRoy} function
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{estbivgeomRoy}}
}
\keyword{distribution}
\keyword{htest}

