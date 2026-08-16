\name{estdweibull}
\alias{estdweibull}
\title{
Estimation of parameters
}
\description{
Estimation of the parameters of the type 1 discrete Weibull distribution
}
\usage{
estdweibull(x, method = "ML", zero = FALSE, eps = 1e-04, nmax=1000)
}
\arguments{
  \item{x}{
the vector of sample values
}
  \item{method}{
\code{"ML"} for the maximum likelihood method; \code{"M"} for the method of moments; \code{"P"} for the method of proportions
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
  \item{eps}{
error threshold for the computation of the moments of the distribution
}
\item{nmax}{
maximum value considered for the numerical computation of the expected value
}
}
\value{
the vector of the estimates of \eqn{q} and \eqn{\beta}
}
\author{
Alessandro Barbiero
}
\seealso{\code{\link{ddweibull}}
}
\examples{
# Ex1
n <- 10
q <- 0.5
beta <- 0.8
x <- rdweibull(n, q, beta)
estdweibull(x, "ML") # maximum likelihood method
# it may return some harmless warnings
# that depend on the optimization function used in the maximization routine
estdweibull(x, "M") # method of moments
estdweibull(x, "P") # method of proportion
# the estimates provided by the three methods may be quite different
# from the true values... and to each other
# change the sample size
n <- 50
q <- 0.5
beta <- 0.8
x <- rdweibull(n, q, beta)
estdweibull(x, "ML") # maximum likelihood method
estdweibull(x, "M") # method of moments
estdweibull(x, "P") # method of proportion
# the estimates should be (on average) closer to the true values
# ...and to each other

# When the estimation methods fail...
# Ex2
# only 1s and 2s
x <- c(1,1,1,1,1,1,2,2,2,2)
estdweibull(x, "ML") # fails!
estdweibull(x, "M") # fails!
estdweibull(x, "P") # fails!

# Ex3
# no 1s
x <- c(2,2,3,4,5,5,5,6,6,8,10)
estdweibull(x, "ML") # works
estdweibull(x, "M") # works
estdweibull(x, "P") # fails!
}
\keyword{htest}
\keyword{distribution}
