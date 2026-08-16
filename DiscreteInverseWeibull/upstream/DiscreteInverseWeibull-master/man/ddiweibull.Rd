\name{Discrete Inverse Weibull}
\alias{ddiweibull}
\alias{pdiweibull}
\alias{qdiweibull}
\alias{rdiweibull}
\title{
The discrete inverse Weibull distribution
}
\description{
Probability mass function, distribution function, quantile function and random generation for the discrete inverse Weibull distribution with parameters \eqn{q} and \eqn{\beta}
}
\usage{
ddiweibull(x, q, beta)
pdiweibull(x, q, beta)
qdiweibull(p, q, beta)
rdiweibull(n, q, beta)
}
\arguments{
  \item{x}{
a vector of quantiles
}
  \item{p}{
a vector of probabilities
}
  \item{q}{
the value of the first parameter, \eqn{q}
}
  \item{beta}{
the value of the second parameter, \eqn{\beta}
}
  \item{n}{
the sample size
}
}
\details{
The discrete inverse Weibull distribution has probability mass function given by \eqn{P(X=x;q,\beta)=q^{(x)^{-\beta}}-q^{(x-1)^{\beta}}}, \eqn{x=1,2,3,\ldots}, \eqn{0<q<1, \beta>0}. Its cumulative distribution function is \eqn{F(x; q, \beta)=q^{x^{-\beta}}}
}
\value{
\code{ddiweibull} gives the probability, \code{pdiweibull} gives the distribution function, \code{qdiweibull} gives the quantile function, and \code{rdiweibull} generates random values. See the reference below for the continuous inverse Weibull distribution.
}
\references{
Dutang, C., Goulet, V., Pigeon, M. (2008) actuar: An R package for actuarial
science, Journal of Statistical Software 25(7): 1-37
}
\examples{
# Ex.1
x<-1:10
q<-0.6
beta<-0.8
ddiweibull(x, q, beta)
t<-qdiweibull(0.99, q, beta)
t
pdiweibull(t, q, beta)
# Ex.2
q<-0.4
beta<-1.7
n<-100
x<-rdiweibull(n, q, beta)
tabulate(x)/sum(tabulate(x))
y<-1:round(max(x))
# compare with
ddiweibull(y, q, beta)
}
\keyword{distribution}

