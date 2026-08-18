\name{FyxbivgeomRoy}
\alias{FyxbivgeomRoy}
\title{
Conditional distribution
}
\description{
Conditional distribution function of \eqn{Y} given \eqn{X=x}
}
\usage{
FyxbivgeomRoy(y, theta1, theta2, theta3, x)
}
\arguments{
  \item{y}{
vector of observations from \eqn{Y}
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
  \item{x}{
value of the conditioning variable \eqn{X}
}
}
\value{
The value of the conditional cumulative distribution function \eqn{F_{Y|x}} in \eqn{y}. Used in \code{\link{rbivgeomRoy}} for conditional sampling
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{EyxbivgeomRoy}}, \code{\link{rbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
# probability that Y<=3 given that X=2:
FyxbivgeomRoy(3, theta1, theta2, theta3, 2)
# the unconditional probability would be
pgeom(3, 1-theta2) # i.e. a geometric distribution with parameter 1-theta2
}
\keyword{distribution}

