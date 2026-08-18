\name{RelbivgeomRoy}
\alias{RelbivgeomRoy}
\title{
Reliability parameter
}
\description{
Stress-strength reliability parameter \eqn{R} for Roy's bivariate geometric model
}
\usage{
RelbivgeomRoy(theta1, theta2, theta3)
}
\arguments{
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
The probability \eqn{R:=P(X\leq Y)} for Roy's bivariate geometric model - see Barbiero (2018) for its computation
}
\references{
Barbiero, A. (2018) Properties and estimation of a bivariate geometric model with locally constant failure rates, submitted

Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{dbivgeomRoy}}, \code{\link{FbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
RelbivgeomRoy(theta1, theta2, theta3)
# theoretical stress-strength reliability parameter R=P(X<=Y)
}
\keyword{distribution}

