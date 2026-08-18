\name{estbivgeomRoy}
\alias{estbivgeomRoy}
\title{
Parameter estimation
}
\description{
Parameter estimation for Roy's bivariate geometric model
}
\usage{
estbivgeomRoy(x, y, method = "LS")
}
\arguments{
  \item{x}{
vector of observations from the first variable \eqn{X}
}
  \item{y}{
vector of observations from the first variable \eqn{y}, same length as \code{x}
}
  \item{method}{
One of the possible estimation methods: "ML" (maximum likelihood), "LS" (least squares), "MMP" (method of moment and poroportion), "M1", "M2", "M3", and "M4" (several variants of the method of moments)
}
}
\value{
a vector of length 3 containing the estimates of \eqn{theta_1}, \eqn{theta_2}, and \eqn{theta_3}
}
\references{
Barbiero, A. (2018) Properties and estimation of a bivariate geometric
model with locally constant failure rates, submitted

Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{dbivgeomRoy}}, \code{\link{minuslogRoy}}
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
hattheta # MLEs
estbivgeomRoy(d[,1], d[,2], "LS")
estbivgeomRoy(d[,1], d[,2], "MMP")
}
\keyword{distribution}
\keyword{htest}

