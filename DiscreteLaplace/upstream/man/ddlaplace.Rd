\name{ddlaplace}
\alias{ddlaplace}
\alias{pdlaplace}
\alias{qdlaplace}
\alias{rdlaplace}
\title{
Probability mass function, cumulative distribution function, quantile function and random generation of the DSL
}
\description{
The function computes the probability mass function, the cumulative distribution function, the quantile function of the DSL and implements random generation.
}
\usage{
ddlaplace(x, p, q)
pdlaplace(x, p, q)
qdlaplace(prob, p, q)
rdlaplace(n, p, q)
}
\arguments{
  \item{x}{
vector of quantiles
}
  \item{p}{
the first parameter \eqn{p} in \eqn{(0,1)} of the SDL
}
  \item{q}{
the second parameter \eqn{q} in \eqn{(0,1)} of the SDL
}
  \item{prob}{
vector of probabilities
}
  \item{n}{
number of observations
}
}
\details{
The pmf of the SDL is given by

\eqn{P(X=x; p, q)=\frac{(1-p)(1-q)}{1-pq}p^x; x=0,1,2,3,\dots}

\eqn{P(X=x; p, q)=\frac{(1-p)(1-q)}{1-pq}q^{|x|}; x=0,-1,-2,-3,\dots}

whereas the cumulative distribution function is given by

\eqn{F(x; p, q)=P(X\leq x)=\frac{(1-p)q^{-\lfloor x\rfloor}}{1-pq},x<0}

\eqn{F(x; p, q)=P(X\leq x)=1-\frac{(1-q)p^{\lfloor x\rfloor+1}}{1-pq},x\geq 0}
}
\value{
\code{ddlaplace} returns the probability of \code{x}; \code{pdlaplace} returns the cumulate probability of \code{x}; \code{qdlaplace} returns the \code{prob}- quantile; \code{rdlaplace} returns a random sample of size \code{n} from DSL.
}
\references{
Tomasz J. Kozubowski, Seidu Inusah (2006) A skew Laplace distribution on integers, \emph{Annals of the Institute of Statistical Mathematics}, 58: 555-571
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{ddlaplace2}}
}
\examples{
# pmf
p<-0.7
q<-0.45
x<--10:10
prob<-ddlaplace(x, p, q)
plot(x, prob, type="h")
prob<-ddlaplace(x, q, p) # swap the parameters
plot(x, prob, type="h")
ddlaplace(-4:4, 1:9/10, 9:1/10) # letting p and q be vectors
# cdf
p<-0.2
q<-0.5
x<-c(-3, -1, pi)
pdlaplace(x, p, q)
# quantile function
p<-0.8
q<-0.4
prob<-c(0.2,0.5,0.8)
x<-qdlaplace(prob, p, q)
x # check
upper<-pdlaplace(x, p, q)
upper
lower<-pdlaplace(x-1, p, q)
lower
lower<=prob & prob<=upper
# random generation
n<-100
p<-0.3
q<-0.5
x<-rdlaplace(n, p, q)
x
t<-table(x)
t
plot(t)
}
\keyword{distribution}

