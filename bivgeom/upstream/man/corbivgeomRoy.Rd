\name{corbivgeomRoy}
\alias{corbivgeomRoy}
\title{
Linear correlation
}
\description{
Linear correlation for Roy's bivariate geometric model
}
\usage{
corbivgeomRoy(theta1, theta2, theta3)
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
the value of Pearson's linear correlation - see Barbiero (2018). The linear correlation for Roy's bivariate geometric distribution is negative (or null, for \eqn{\theta_3=1}) for any feasible choice of its parameters
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.

Barbiero, A. (2018) Properties and estimation of a bivariate geometric model with locally constant failure rates, submitted
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{dbivgeomRoy}}
}
\examples{
corbivgeomRoy(0.3,0.7,0.5)
}
\keyword{distribution}

