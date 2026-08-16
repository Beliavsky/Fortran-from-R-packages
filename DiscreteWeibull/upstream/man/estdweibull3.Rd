\name{estdweibull3}
\alias{estdweibull3}
\title{
Estimation of parameters
}
\description{
Estimation of the parameters of the type 3 discrete Weibull distribution
}
\usage{
estdweibull3(x, method = "P", eps = 1e-04)
}
\arguments{
  \item{x}{
the vector of sample values
}
  \item{method}{
\code{"ML"} for the maximum likelihood method; \code{"M"} for the method of moments; \code{"P"} for the method of proportions
}
  \item{eps}{
error threshold for the computation of the moments of the distribution
}
}
\value{
the vector of the estimates of \eqn{c} and \eqn{\beta}
}
\author{
Alessandro Barbiero
}
\seealso{\code{\link{ddweibull3}}
}
\examples{
# Ex1
x <- c(0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,3,3,4,6)
estdweibull3(x, "P")
estdweibull3(x, "ML")
estdweibull3(x, "M")
# Ex 2
n <- 20
c <- 1/3
beta <- 2/3
x <- rdweibull3(n, c, beta)
estdweibull3(x, "P")
par <- estdweibull3(x, "ML")
par
-loglikedw3(par, x)
par <- estdweibull3(x, "M")
par
lossdw3(par, x)
n <- 50
x <- rdweibull3(n, c, beta)
estdweibull3(x, "P")
estdweibull3(x, "ML")
estdweibull3(x, "M")
n <- 100
x <- rdweibull3(n, c, beta)
estdweibull3(x, "P")
estdweibull3(x, "ML")
estdweibull3(x, "M")
# Ex 3: a piece of simulation study
nSim <- 50
n <- 50
c <- 0.2
beta <- 0.7
par <- matrix(0, nSim, 2)
for(i in 1:nSim)
{
x <- rdweibull3(n, c, beta)
par[i,] <- estdweibull3(x, "ML")
}
op <- par(mfrow = c(1,2))
boxplot(par[,1], xlab=expression(hat(c)[ML]))
abline(h = c)
boxplot(par[,2], xlab=expression(hat(beta)[ML]))
abline(h = beta)
op <- par()
}
\keyword{htest}
\keyword{distribution}
