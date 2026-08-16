\name{loglikediw}
\alias{loglikediw}
\title{
likelihood function
}
\description{
Log-likelihood function of the discrete inverse Weibull
}
\usage{
loglikediw(x, q, beta)
}
\arguments{
  \item{x}{
a vector of sample values
}
  \item{q}{
the value of the \eqn{q} parameter
}
  \item{beta}{
the value of the \eqn{\beta} parameter
}
}
\value{
the value of the log-likelihood function (changed in sign) of the discrete inverse Weibull distribution with parameters \eqn{q} and \eqn{\beta} computed on a sample x
}
\seealso{
\code{\link{heuristic}}
}
\examples{
n<-100
q<-0.4
beta<-2
x<-rdiweibull(n, q, beta)
# loglikelihood function (changed in sign) computed on the true values
loglikediw(x, q, beta)
par<-estdiweibull(x, method="H")
par
# loglikelihood function (changed in sign) computed on the ML estimates
loglikediw(x, par[1], par[2])
# it should be smaller than before...
}
\keyword{distribution, htest}
