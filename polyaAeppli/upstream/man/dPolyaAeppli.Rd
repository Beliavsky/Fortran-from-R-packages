\name{PolyaAeppli}
\alias{dPolyaAeppli}
\alias{pPolyaAeppli}
\alias{qPolyaAeppli}
\alias{rPolyaAeppli}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Polya-Aeppli 
}
\description{
Density, distribution function, quantile function and random generation for the Polya-Aeppli distribution with parameters lambda and prob. 
}
\usage{
dPolyaAeppli(x, lambda, prob, log = FALSE)
pPolyaAeppli(q, lambda, prob, lower.tail = TRUE, log.p = FALSE)
qPolyaAeppli(p, lambda, prob, lower.tail = TRUE, log.p = FALSE)
rPolyaAeppli(n, lambda, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles
}
  \item{q}{
vector of quantiles
}
  \item{p}{
vector of probabilities
}
  \item{n}{
number of random variables to return
}
  \item{lambda}{
a vector of non-negative Poisson parameters
}
  \item{prob}{
a vector of geometric parameters between 0 and 1
}
  \item{log, log.p}{
logical; if TRUE, probabilities p are given as log(p)
}
  \item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}
}
}
\details{
A Polya-Aeppli, or geometric compound Poisson, random variable is the sum of a Poisson number of identically and independently distributed shifted geometric random variables.  Its distribution (with \code{lambda}\eqn{= \lambda}, \code{prob}\eqn{= p}) has density 
\deqn{Prob(X = x) = e^(-\lambda)} for \eqn{x = 0};
\deqn{Prob(X = x) = e^(-\lambda) \sum_{n = 1}^y (\lambda^n)/(n!) choose(y - 1, n - 1) p^(y - n) (1 - p)^n} 
for \eqn{x = 1, 2, \ldots}.  

If an element of x is not integer, the result of \code{dPolyaAeppli} is zero, with a warning. 

The quantile is right continuous: \code{qPolyaAeppli(p, lambda, prob)} is the smallest integer \eqn{x} such that \eqn{P(X \le x) \ge p}.

Setting \code{lower.tail = FALSE} enables much more precise results when the default, lower.tail = TRUE would return 1, see the example below.

}
\value{
\code{dPolyaAeppli} gives the (log) density, \code{pPolyaAepploi} gives the (log) distribution function, \code{qPolyaAeppli} gives the quantile function, and \code{rPolyaAeppli} generates random deviates.

Invalid \code{lambda} or \code{prob} will terminate with an error message.

}
\references{
Johnson NL, Kotz S, Kemp AW (1992). \emph{Univariate Discrete Distributions.} 2nd edition. Wiley, New York.

Nuel G (2008). \emph{Cumulative distribution function of a geometeric Poisson distribution.} Journal of Statistical Computation and Simulation, \bold{78}(3), 385-394.}
\author{
Conrad Burden
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
