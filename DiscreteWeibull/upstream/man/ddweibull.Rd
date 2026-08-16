\name{Discrete Weibull (Type 1)}
\alias{ddweibull}
\alias{pdweibull}
\alias{qdweibull}
\alias{rdweibull}
\title{
The type 1 discrete Weibull distribution
}
\description{
Probability mass function, distribution function, quantile function and random generation for the discrete Weibull distribution with parameters \eqn{q} and \eqn{\beta}
}
\usage{
ddweibull(x, q, beta, zero = FALSE)
pdweibull(x, q, beta, zero = FALSE)
qdweibull(p, q, beta, zero = FALSE)
rdweibull(n, q, beta, zero = FALSE)
}
\arguments{
  \item{x}{
vector of quantiles
}
  \item{p}{
vector of probabilities
}
  \item{q}{
first parameter
}
  \item{beta}{
second parameter
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
  \item{n}{
sample size
}
}
\details{
The discrete Weibull distribution has probability mass function given by \eqn{P(X=x;q,\beta)=q^{(x-1)^{\beta}}-q^{x^{\beta}}}, \eqn{x=1,2,3,\ldots}, if \code{zero}=\code{FALSE}; or \eqn{P(X=x;q,\beta)=q^{x^{\beta}}-q^{(x+1)^{\beta}}}, \eqn{x=0,1,2,\ldots}, if \code{zero}=\code{TRUE}. The cumulative distribution function is \eqn{F(x;q,\beta)=1-q^{x^{\beta}}} if \code{zero}=\code{FALSE}; \eqn{F(x;q,\beta)=1-q^{(x+1)^{\beta}}} otherwise 
}
\value{
\code{ddweibull} gives the probability function, \code{pdweibull} gives the distribution function, \code{qdweibull} gives the quantile function, and \code{rdweibull} generates random values.
}
\author{
Alessandro Barbiero
}
\examples{
# Ex.1
x <- 1:10
q <- 0.6
beta <- 0.8
ddweibull(x, q, beta)
t <- qdweibull(0.99, q, beta)
t
pdweibull(t, q, beta)
# 
x <- 0:10
ddweibull(x, q, beta, zero=TRUE)
t <- qdweibull(0.99, q, beta, zero=TRUE)
t
pdweibull(t, q, beta, zero=TRUE)

# Ex.2
q <- 0.4
beta <- 0.7
n <- 100
x <- rdweibull(n, q, beta)
tabulate(x)/sum(tabulate(x))
y <- 1:round(max(x))
# compare with
ddweibull(y, q, beta)
}
\keyword{distribution}

