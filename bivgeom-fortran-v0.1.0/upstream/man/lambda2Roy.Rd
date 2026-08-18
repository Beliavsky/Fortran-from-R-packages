\name{lambda2Roy}
\alias{lambda2Roy}
\title{
Bivariate failure rate
}
\description{
Bivariate failure rate \eqn{\lambda_2}
}
\usage{
lambda2Roy(x, y, theta1, theta2, theta3)
}
\arguments{
  \item{x}{
observation from the first variable
}
  \item{y}{
observation from the second variable
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
Value of the bivariate failure rate \eqn{\lambda_2} for Roy's bivariate geometric model (Roy, 1993)
}
\details{
It is defined as \eqn{P(X\geq x,Y=y)/P(X\geq x,Y\geq y)}. For this model, \eqn{\lambda_2(x,y)=1-\theta_2\theta_3^x}
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{lambda1Roy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
# bivariate failure rate lambda 2
# computed in x=1, y=2
x <- 1
y <- 2
lambda2Roy(x,y,theta1,theta2,theta3)
}
\keyword{distribution}

