\name{gb2}
\Rdversion{2.1}
\alias{gb2}
\alias{dgb2}
\alias{pgb2}
\alias{qgb2}
\alias{rgb2}
\title{The Generalized Beta Distribution of the Second Kind}
\description{Density, distribution function, quantile function and random generation for the Generalized beta distribution of the second kind with parameters \code{a}, \code{b}, \code{p} and \code{q}. 
}
\usage{
dgb2(x, shape1, scale, shape2, shape3)
pgb2(x, shape1, scale, shape2, shape3)
qgb2(prob, shape1, scale, shape2, shape3)
rgb2(n, shape1, scale, shape2, shape3)
}
\arguments{
  \item{x}{numeric; vector of quantiles.}
  \item{shape1}{numeric; positive parameter.}
  \item{scale}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
  \item{prob}{numeric; vector of probabilities.}
  \item{n}{numeric; number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
}
\details{
The Generalized Beta distribution of the second kind with parameters \code{shape1} \eqn{= a}, \code{scale} \eqn{= b},  \code{shape2} \eqn{= p} and \code{shape3} \eqn{= q} has density
  \deqn{f(x)=\frac{a(x/b)^{ap-1}}{bB(p,q)(1+(x/b)^{a})^{p+q}}}{a(x/b)^(ap-1)/bB(p,q)(1+(x/b)^(a))^(p+q), x \ge 0}
  for \eqn{a > 0}, \eqn{b > 0}, \eqn{p > 0} and \eqn{q > 0}, where  \eqn{B(p,q)} is the Beta function (\code{\link{beta}}). If \code{Z} follows a Beta distribution with parameters \eqn{p} and \eqn{q} and 
   \deqn{y = \frac{z}{1-z},} then \deqn{x = b * y^{1/a}} follows the GB2 distribution.
}
\value{
\code{dgb2} gives the density, \code{pgb2} the distribution
  function, \code{qgb2} the quantile function, and \code{rgb2}
  generates random deviates.
}
\references{
 Kleiber, C. and Kotz, S. (2003) 
 \emph{Statistical Size Distributions in Economics and Actuarial Sciences}, chapter 6.
 Wiley, Ney York.
 
 McDonald, J. B. (1984)
 Some generalized functions for the size distribution of income. 
 \emph{Econometrica}, \bold{52}, 647--663.
}
\author{Monique Graf}
\seealso{ \code{\link{beta}} for the Beta function and  \code{\link{dbeta}} for the Beta distribution.
}
\examples{
a <- 3.9
b <- 18873
p <- 0.97
q <- 1.03
x <- qgb2(0.6, a, b, p, q)
y <- dgb2(x, a, b, p, q)
}
\keyword{distribution}

