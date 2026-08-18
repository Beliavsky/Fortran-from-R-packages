\name{loglikgeomRoy}
\alias{loglikgeomRoy}
\title{
Log-likelihood function
}
\description{
Negative log-likelihood function for Roy's bivariate geometric model
}
\usage{
loglikgeomRoy(par, x, y)
}
\arguments{
  \item{par}{
a vector containing the values of the three parameters \eqn{\theta_1}, \eqn{\theta_2}, and \eqn{\theta_3}
}
  \item{x}{
numeric vector of sample \eqn{x}-values (non-negative integers)
}
  \item{y}{
numeric vector of sample \eqn{x}-values (non-negative integers), same length as \code{x}
}
}
\value{
Value of the negative log-likelihood function
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{dbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
# random sample of size n=1000:
set.seed(12345)
n <- 1000
d <- rbivgeomRoy(n, theta1, theta2, theta3)
# parameter estimation, using the different proposed methods:
hattheta <- estbivgeomRoy(d[,1], d[,2], "ML")
loglikgeomRoy(hattheta, x=d[,1], y=d[,2])
# negative value of the (maximized) log-likelihood function
}
\keyword{distribution}
\keyword{htest}

