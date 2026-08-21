\name{Moments}
\Rdversion{2.1}
\alias{Moments}
\alias{moment.gb2}
\alias{incompl.gb2}
\alias{el.gb2}
\alias{vl.gb2}
\alias{sl.gb2}
\alias{kl.gb2}
\title{
Moments and Other Properties of a GB2 Random Variable
}
\description{
These functions calculate the moments of order \code{k} and incomplete moments of order \code{k} of a GB2 random variable \eqn{X} as well as the expectation,
the variance, the kurtosis and the skewness of \eqn{log(X)}.
}
\usage{
moment.gb2(k, shape1, scale, shape2, shape3)
incompl.gb2(x, k, shape1, scale, shape2, shape3)
el.gb2(shape1, scale, shape2, shape3)
vl.gb2(shape1, shape2, shape3)
sl.gb2(shape2, shape3)
kl.gb2(shape2, shape3)
}
\arguments{
  \item{x}{numeric; vector of quantiles.}
  \item{k}{numeric; order of the moment.}
  \item{shape1}{numeric; positive parameter.}
  \item{scale}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
}
\details{
Let \eqn{X} be a random variable following a GB2 distribution with parameters \code{shape1} \eqn{= a}, \code{scale} \eqn{= b},  \code{shape2} \eqn{= p} and \code{shape3} \eqn{= q}. 
Moments and incomplete moments of \eqn{X} exist only for \eqn{-ap \le k \le aq}{-ap <= k <= aq}. Moments are given by
\deqn{E(X^k) = {b}^{k} \frac{\Gamma (p+k/a) \Gamma (q-k/a)}{\Gamma (p) \Gamma (q)}}{E(X^k) = b^k Gamma(p+k/a) Gamma(q-k/a) / Gamma(p) Gamma(q).} 

This expression, when considered a function of \code{k}, can be viewed as the moment-generating function of \eqn{Y=log(X)}. Thus, it is useful to compute the moments of \eqn{log(X)}, 
which are needed for deriving, for instance, the Fisher information matrix of the GB2 distribution. Moments of \eqn{log(X)} exist for all \code{k}.
}
\value{
\code{moment.gb2} gives the moment of order \code{k}, 
\code{incompl.gb2} gives the incomplete moment of order \code{k},
\code{El.gb2} gives the expectation of \eqn{log(X)}, 
\code{vl.gb2} gives the variance of \eqn{log(X)}, 
\code{sl.gb2} gives the skewness of \eqn{log(X)}, 
\code{kl.gb2} gives the kurtosis of \eqn{log(X)}.
}
\references{
 Kleiber, C. and Kotz, S. (2003) 
 \emph{Statistical Size Distributions in Economics and Actuarial Sciences}, chapter 6.
 Wiley, Ney York.
}
\author{
Monique Graf
}
\seealso{ 
\code{\link{gamma}} for the Gamma function and related functions (\code{digamma}, \code{trigamma} and \code{psigamma}).
}
\examples{
a <- 3.9
b <- 18873
p <- 0.97
q <- 1.03
k <- 2
x <- qgb2(0.6, a, b, p, q)
moment.gb2(k, a, b, p, q)
incompl.gb2(x, k, a, b, p, q)
vl.gb2(a, p, q)
kl.gb2(p, q)
}
\keyword{distribution}

