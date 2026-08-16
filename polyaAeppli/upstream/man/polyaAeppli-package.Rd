\name{polyaAeppli-package}
\alias{polyaAeppli-package}
\alias{polyaAeppli}
\docType{package}
\title{
Implementation of the Polya-Aeppli Distribution
}
\description{
Functions for evaluating the mass density, cumulative distribution function, quantile function and random variate generation for the Polya-Aeppli distribution, also known as the geometric compound Poisson distribution. \cr
More information on the implementation of \bold{polyaAeppli} can be found at Conrad J. Burden (2014) <arXiv:1406.2780>.
}

\details{
\tabular{ll}{
Package: \tab polyaAeppli\cr
Type: \tab Package\cr
Version: \tab 2.0.2\cr
Depends: \tab R (>= 3.0.0)\cr
Date: \tab 2020-04-21\cr
License: \tab GPL(>=2)\cr
}
Consistent with the conventions used in R package stats, this implementation of the Polya-Aeppli distribution comprises the four functions \cr

\code{dPolyaAeppli(x, lambda, prob, log = FALSE)}\cr
\code{pPolyaAeppli(q, lambda, prob, lower.tail = TRUE, log.p = FALSE)}\cr
\code{qPolyaAeppli(p, lambda, prob, lower.tail = TRUE, log.p = FALSE)} \cr
\code{rPolyaAeppli(n, lambda, prob)}
}

\author{
Conrad Burden

Maintainer: conrad.burden@anu.edu.au
}
\references{
Johnson NL, Kotz S, Kemp AW (1992). \emph{Univariate Discrete Distributions.} 2nd edition. Wiley, New York.

Nuel G (2008). \emph{Cumulative distribution function of a geometeric Poisson distribution.} Journal of Statistical Computation and Simulation, \bold{78}(3), 385-394.
}
\examples{
lambda <- 8
prob <- 0.2
## Plot histogram of random sample
PAsample <- rPolyaAeppli(10000, lambda, prob)  
maxPA <- max(PAsample)
hist(PAsample, breaks=(0:(maxPA + 1)) - 0.5, freq=FALSE, 
  		xlab = "x", ylab = expression(P[X](x)), main="", border="blue")
## Add plot of density function
x <- 0:maxPA
points(x, dPolyaAeppli(x, lambda, prob), type="h", lwd=2)

lambda <- 4000
prob <- 0.005
qq <- 0:10000
## Plot log of the extreme lower tail p-value
log.pp <- pPolyaAeppli(qq, lambda, prob, log.p=TRUE)
plot(qq, log.pp, type = "l", ylim=c(-lambda,0), 
  				xlab = "x", ylab = expression("log Pr(X " <= "x)"))
## Plot log of the extreme upper tail p-value
log.1minuspp <- pPolyaAeppli(qq, lambda, prob, log.p=TRUE, lower.tail=FALSE)
points(qq, log.1minuspp, type = "l", col = "red")
legend("topright", c("lower tail", "upper tail"), 
				col=c("black", "red"), lty=1, bg="white") 
}