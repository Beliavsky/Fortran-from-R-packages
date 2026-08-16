\name{ddlaplace2}
\alias{ddlaplace2}
\alias{pdlaplace2}
\alias{qdlaplace2}
\alias{rdlaplace2}
\alias{palaplace2}
\title{
Probability mass function of the ADSL
}
\description{
The function computes the probability mass function, the cumulative distribution function, the quantile function of the ADSL and provides random generation of samples from the same model
}
\usage{
ddlaplace2(x, p, q)
palaplace2(x, p, q)
pdlaplace2(x, p, q)
qdlaplace2(prob, p, q)
rdlaplace2(n, p, q)
}
\arguments{
  \item{x}{
vector of quantiles
}
  \item{p}{
the first parameter \eqn{p} in \eqn{(0,1)} of the ADSL
}
  \item{q}{
the second parameter \eqn{q} in \eqn{(0,1)} of the ADSL
}
  \item{prob}{
vector of probabilities
}
  \item{n}{
number of observations
}
}
\details{
The probability mass funtion of the ADSL distribution is given by:

\eqn{P(X=x;p,q)=\frac{\log p}{\log (pq)}q^{-(x+1)}(1-q)} for \eqn{x=\dots, -2, -1}

and

\eqn{P(X=x;p,q)=\frac{\log q}{\log (pq)}p^{x}(1-p)} for \eqn{x=0, 1, 2, \dots}

Its cumulative distribution function is:

\eqn{F(x;p,q)=\frac{\log p}{\log (pq)}q^{-(\lfloor x \rfloor+1)}} for \eqn{x<0}

and

\eqn{F(x;p,q)=1-\frac{\log q}{\log (pq)}p^{(\lfloor x \rfloor+1)}} for \eqn{x\geq 0}
}
\value{
\code{ddlaplace2} returns the probability of \code{x}; \code{pdlaplace2} returns the cumulate probability of \code{x}; \code{qdlaplace2} returns the \code{prob}- quantile; \code{rdlaplace2} returns a random sample of size \code{n} from ADSL.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{ddlaplace}}
}
\examples{
# pmf
p <- 0.7
q <- 0.45
x <- -10:10
prob <- ddlaplace2(x, p, q)
plot(x, prob, type="h")
# swap the parameters
prob <- ddlaplace2(x, q, p)
plot(x, prob, type="h")
# letting p and q be vectors...
ddlaplace2(-4:4, 1:9/10, 9:1/10)
# cdf
pdlaplace2(x, p, q)
pdlaplace2(pi, p, q)
pdlaplace2(floor(pi), p, q)
# quantile function
qdlaplace(1:9/10, p, q)
# random generation
y <- rdlaplace2(n=1000, p, q)
plot(table(y))
}
\keyword{distribution}

