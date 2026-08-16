\name{Discrete Weibull (Type 3)}
\alias{ddweibull3}
\alias{pdweibull3}
\alias{qdweibull3}
\alias{rdweibull3}
\alias{hdweibull3}
\title{
The type 3 discrete Weibull distribution
}
\description{
Probability mass function, distribution function, quantile function, random generation, and hazard function for the type 3 discrete Weibull distribution with parameters \eqn{c} and \eqn{\beta}
}
\usage{
ddweibull3(x, c, beta)
pdweibull3(x, c, beta)
qdweibull3(p, c, beta)
rdweibull3(n, c, beta)
hdweibull3(x, c, beta)
}
\arguments{
  \item{x}{
vector of values/quantiles
}
  \item{p}{
vector of probabilities
}
  \item{c}{
first parameter
}
  \item{beta}{
second parameter
}
  \item{n}{
sample size
}
}
\details{
The type 3 discrete Weibull distribution is characterized by the following cumulative distribution function: \eqn{F(x;c,\beta)=1-\exp(-c\sum_{j=0}^{x+1}j^{\beta})}, for \eqn{x=0,1,2,\dots}, with \eqn{c>0} and \eqn{\beta>=-1}. 
}
\value{
\code{ddweibull3} gives the probability function, \code{pdweibull3} gives the distribution function, \code{qdweibull3} gives the quantile function, \code{hdweibull3} gives the hazard function, and \code{rdweibull} generates random values.
}
\author{
Alessandro Barbiero
}
\examples{
# ddweibull3
x <- 0:10
c <- 0.3
beta <- 0.75
p <- ddweibull3(x, c, beta)
p
plot(x, p, type="b", ylab=expression(P(X)==x))
# pdweibull3
x <- 0:10
c <- 0.5
beta <- 0.5
p <- pdweibull3(x, c, beta)
p
cumsum(ddweibull3(x, c, beta))
plot(x, p, type="s", ylab=expression(P(X<=x)))
# qdweibull3
p <- c(1:9)/10
p
c <- 0.1
beta <- 0.5
qdweibull3(p, c, beta)
pdweibull3(10, c, beta)
pdweibull3(9, c, beta)
# rdweibull3
n <- 20
c <- 0.25
beta <- -0.25
x <- rdweibull3(n, c, beta)
x
beta <- 0
x <- rdweibull3(n, c, beta)
x
beta <- 0.25
x <- rdweibull3(n, c, beta)
x
n <- 1000
x <- rdweibull3(n, c, beta)
obs <- c(sum(x==0), tabulate(x))
obs <- obs/sum(obs)
theo <- ddweibull3(min(x):max(x), c, beta)
barplot(rbind(obs, theo), beside=TRUE, names.arg=min(x):max(x),
ylab="relative frequency/probability", col=1:2)
legend(24, 0.1, c("observed", "theoretical"), pch=15, col=1:2)
#hdweibull3
x <- 0:15
c <- 0.5
hn<-hdweibull3(x, c, beta = -0.5)
h0<-hdweibull3(x, c, beta = 0)
hp<-hdweibull3(x, c, beta = 0.5)
plot(x, hn, type="b", ylim = c(0, 1), ylab="hazard rate")
points(x, h0, type = "b", col=2)
points(x, hp, type = "b", col=3)
legend(11, 0.5, c("beta<0", "beta=0", "beta>0"), col=1:3, pch=21)
}
\keyword{distribution}

