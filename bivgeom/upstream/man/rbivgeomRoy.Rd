\name{rbivgeomRoy}
\alias{rbivgeomRoy}
\title{
Pseudo-random generation
}
\description{
Generation of pseudo-random values from Roy's bivariate geometric model
}
\usage{
rbivgeomRoy(n, theta1, theta2, theta3)
}
\arguments{
  \item{n}{
a positive integer, corresponding to the sample size
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
A \eqn{n\times 2} numeric matrix containing the bivariate sample values 
}
\references{
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
# random sample of size n=1000:
set.seed(12345)
n <- 1000
d <- rbivgeomRoy(n, theta1, theta2, theta3)
# joint frequency distribution:
table(d[,1],d[,2])
}
\keyword{distribution}
\keyword{datagen}

