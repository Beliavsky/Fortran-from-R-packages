\name{Edweibull3}
\alias{Edweibull3}
\alias{E2dweibull3}
\title{
Expected values
}
\description{
First and second order moments for the type 3 discrete Weibull distribution
}
\usage{
Edweibull3(c, beta, eps = 1e-04)
E2dweibull3(c, beta, eps = 1e-04)
}
\arguments{
  \item{c}{
first parameter
}
  \item{beta}{
second parameter
}
  \item{eps}{
error threshold for the numerical computation of the expected value
}
}
\details{
The expected values are numerically computed considering a truncated support: integer values smaller than or equal to
\eqn{2F^{-1}(1-eps;c,\beta))}, where \eqn{F^{-1}} is the inverse of the cumulative distribution function (implemented by the function \code{\link{qdweibull3}})
}
\value{
the (approximate) expected values of the discrete Weibull distribution:
\code{Edweibull3} gives the first order moment, \code{E2dweibull3} the second order moment
}
\author{
Alessandro Barbiero
}
\examples{
c <- 0.4
beta <- 0.25
Edweibull3(c,beta)
c <- 0.4
beta <- -0.75
Edweibull3(c, beta) # may require too much time
Edweibull3(c, beta, eps=0.001) # try with a smaller eps->worse approximation
c <- rep(0.1, 11)
beta <- (0:10)/10
Edweibull3(c, beta)
c <- rep(0.5, 11)
beta <- (-5:5)/10
Edweibull3(c,beta)
# E2dweibull3
c <- 0.4
beta <- 0.25
E2dweibull3(c, beta)
c <- rep(0.1, 11)
beta <- (0:10)/10
Edweibull3(c, beta)
c <- rep(0.8, 11)
beta <- (-5:5)/11
E2dweibull3(c, beta)
}
\keyword{distribution}