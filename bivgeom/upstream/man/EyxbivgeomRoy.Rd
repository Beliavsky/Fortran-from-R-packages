\name{EyxbivgeomRoy}
\alias{EyxbivgeomRoy}
\title{
Conditional moment
}
\description{
Conditional moment of \eqn{Y} given \eqn{X=x} for Roy's bivariate geomtric model
}
\usage{
EyxbivgeomRoy(theta1, theta2, theta3, x)
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
  \item{x}{
value of the conditioning variable \eqn{X}
}
}
\value{
Value of the conditional moment of \eqn{Y} given \eqn{X=x}
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}

\seealso{
\code{\link{FyxbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
EyxbivgeomRoy(theta1, theta2, theta3, 2)
}
\keyword{distribution}

