\name{Edweibull}
\alias{Edweibull}
\alias{E2dweibull}
\alias{Vdweibull}
\alias{ERdweibull}
\title{
Expected values
}
\description{
First and second order moments, variance and expected value of the reciprocal for the type 1 discrete Weibull distribution
}
\usage{
Edweibull(q, beta, eps = 1e-04, nmax = 1000, zero = FALSE)
E2dweibull(q, beta, eps = 1e-04, nmax = 1000, zero = FALSE)
Vdweibull(q, beta, eps = 1e-04, nmax = 1000, zero = FALSE)
ERdweibull(q, beta, eps = 1e-04, nmax = 1000)
}
\arguments{
  \item{q}{
first parameter
}
  \item{beta}{
second parameter
}
  \item{eps}{
error threshold for the numerical computation of the expected value
}
  \item{nmax}{
maximum value considered for the numerical approximate computation of the expected value;
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
}
\details{
The expected value is numerically computed considering a truncated support: integer values smaller than or equal to \eqn{2F^{-1}(1-eps;q,\beta)} are considered, where \eqn{F^{-1}} is the inverse of the cumulative distribution function (implemented by the function \code{\link{qdweibull}}). However, if such value is greater than \code{nmax}, the expected value is computed recalling the formula of the expected value of the corresponding continuous Weibull distribution (see the reference), adding 0.5. Similar arguments apply to the other moments.
}
\value{
the (approximate) expected values of the discrete Weibull distribution:
\code{Edweibull} gives the first order moment, \code{E2dweibull} the second order moment, \code{Vdweibull} the variance, \code{ERdweibull} the expected value of the reciprocal (only if \code{zero} is \code{FALSE})
}
\author{
Alessandro Barbiero
}
\examples{
q <- 0.75
beta <- 1.25
Edweibull(q, beta)
E2dweibull(q, beta)
Vdweibull(q, beta)
ERdweibull(q, beta)
# if beta=0.75...
beta <- 0.75
Edweibull(q, beta)
Edweibull(q, beta, nmax=100)
# here above, the approximation through the continuous model intervenes
# if beta=1...
beta <- 1
Edweibull(q, beta)
# which equals...
1/(1-q)
}
\references{
M. S. A. Khan, A. Khalique, and A. M. Abouammoh (1989) On estimating parameters in a discrete Weibull distribution, IEEE Transactions on Reliability, 38(3), pp. 348-350
}
\keyword{distribution}