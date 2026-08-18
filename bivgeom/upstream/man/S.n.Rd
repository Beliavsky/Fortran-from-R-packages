\name{S.n}
\alias{S.n}
\title{
Empirical joint survival function
}
\description{
Empirical joint survival function
}
\usage{
S.n(x, X)
}
\arguments{
  \item{x}{
matrix with two columns of non-negative integer values where the empirical joint survival function is computed
}
  \item{X}{
matrix with two columns corresponding to the full observed sample
}
}
\value{
value of the empirical joint survival function \eqn{\hat{S}_X(x)}
}
\references{
Roy, D. (1993) Reliability measures in the discrete bivariate set-up and related characterization results for a bivariate geometric distribution, Journal of Multivariate Analysis 46(2), 362-373.
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{estbivgeomRoy}}
}
\examples{
theta1 <- 0.5
theta2 <- 0.7
theta3 <- 0.9
set.seed(12345)
n <- 1000
d <- rbivgeomRoy(n, theta1, theta2, theta3)
S.n(cbind(1,1),d) # empirical sf
# compare it with the theoretical
SbivgeomRoy(1,1,theta1,theta2,theta3)
}
\keyword{distribution}
\keyword{nonparametric}


