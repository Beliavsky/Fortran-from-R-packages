\name{FbivgeomRoy}
\alias{FbivgeomRoy}
\title{
Joint distribution function
}
\description{
Joint cumulative distribution function for Roy's bivariate geometric model
}
\usage{
FbivgeomRoy(x, y, theta1, theta2, theta3)
}
\arguments{
  \item{x}{
vector of values for the first variable \eqn{X}
}
  \item{y}{
vector of values for the second variable \eqn{Y} 
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
The probability \eqn{F(x,y):=P(X\leq x,Y\leq x)}
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{dbivgeomRoy}}, \code{\link{SbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
# probability that X<=2 and Y<=3:
FbivgeomRoy(2, 3, theta1, theta2, theta3)
}
\keyword{distribution}

