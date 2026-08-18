\name{dbivgeomRoy}
\alias{dbivgeomRoy}
\title{
Joint probability mass function
}
\description{
Joint probability mass function for Roy's bivariate geometric model
}
\usage{
dbivgeomRoy(x, y, theta1, theta2, theta3)
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
Value of the probability \eqn{p(x,y):=P(X=x,Y=y)}. 
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
dbivgeomRoy(x=2, y=0, theta1=0.7, theta2=0.2, theta3=0.8)
dbivgeomRoy(0:5, y=0, theta1=0.7, theta2=0.2, theta3=0.8)
# these are p(0,0), p(1,0), ..., p(5,0)
dbivgeomRoy(0:2, 1:3, theta1=0.7, theta2=0.2, theta3=0.8)
# these are p(0,1), p(1,2), p(2,3)
}
\keyword{distribution}

