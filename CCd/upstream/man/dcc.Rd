\name{dcc}
\alias{dcc}
\alias{pcc}
\alias{qcc}
\title{
Density, distribution function and quantile function of the CC distribution
}

\description{
Density, distribution function and quantile function of the CC distribution.
}

\usage{
dcc(y, mu = 0, lambda, logged = FALSE)
pcc(y, mu = 0, lambda)
qcc(p, mu, lambda)
}

\arguments{
\item{y}{
A vector with integer values.
}
\item{p}{
A vector with probabilities.
}
\item{mu}{
The value of the location parameter \eqn{\mu}.
}
\item{lambda}{
The value of the scale parameter \eqn{\lambda}.
}
\item{logged}{
Should the logarithm of the density be returned (TRUE) or not (FALSE)?
}
}

\details{
The density of the CC distribution is computed. The probability mass function of the CC distribution (Papadatos, 2022) is given by
\eqn{P(X=k)=\dfrac{\tanh{(\lambda \pi)}}{\pi}\dfrac{\lambda}{\lambda^2+\kappa^2}.}

The cumulative distribution function of the CC distribution is computed. We explore the property of the CC distribution that \eqn{P(X=-\kappa)=P(X=\kappa)}, where \eqn{\kappa>0}, to compute the cumulative distribution.

As for the quantile function we use the \code{\link{optimize}} function to find the integer whose cumulative probability matches the given probability. So, basically, the qcc() works with left tailed probabilities.
}

\value{
\bold{dcc} returns a vector with the (logged) density values, the (logged) probabilities for each value of y., \bold{pcc} returns a vector with the cumulative probabilities, while \bold{qcc} returns a vector with integer numbers.
}

\references{
Papadatos N. (2022). The characteristic function of the discrete Cauchy distribution In Memory of T. Cacoullos. Journal of Statistical Theory and Practice, 16(3): 47.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{dcc}, \link{cc.mle} }
}

\examples{
x <- round( rcauchy(100, 3, 10) )
mod <- cc.mle(x)
y <- dcc(x, mod$param[1], mod$param[3])

pcc(x[1:5], mod$param[1], mod$param[3])
}

