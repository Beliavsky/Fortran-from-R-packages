\name{SbivgeomRoy}
\alias{SbivgeomRoy}
\title{
Joint survival function
}
\description{
Joint survival function for Roy's bivariate geometric model
}
\usage{
SbivgeomRoy(x, y, theta1, theta2, theta3)
}
\arguments{
  \item{x}{
vector of observations from the first variable \eqn{X}
}
  \item{y}{
vector of observations from the second variable \eqn{Y} (same length as \code{x})
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
The probability \eqn{P(X\geq x, Y\geq y}. For this model it is equal to \eqn{S(x,y)=\theta_1^x\theta_2^y\theta_3^{xy}}
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{FbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
# probability that X>=2 and Y>=3:
SbivgeomRoy(2, 3, theta1, theta2, theta3)
}
\keyword{distribution}

