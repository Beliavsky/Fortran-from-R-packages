\name{Generalized Binomial}
\alias{Generalized Binomial}
\alias{dgbinom}
\alias{pgbinom}
\alias{qgbinom}
\alias{rgbinom}
\title{
The Generalized Binomial Distribution
}
\description{
Density, distribution function, quantile function and random generation for the generalized binomial distribution with parameter vectors \code{size} and \code{prob}.
}
\usage{
dgbinom(x, size, prob, log = FALSE)
pgbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
qgbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)
rgbinom(N, size, prob)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
\item{p}{
vector of probabilities.
}
\item{N}{
number of observations.
}
  \item{size}{
vector of the number of trials for each type.
}
  \item{prob}{
vector of the success probabilities  for each type.
}
  \item{log,log.p}{
logical; if TRUE probabilities p are given as log(p).
}
\item{lower.tail}{
logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}{P[X\le x]}, otherwise, \eqn{P[X>x]} .
}
}
\details{
The generalized binomial distribution with \code{size}=\eqn{c(n_{1},\dots ,n_{r})}{c(n1,\dots ,nr)} and \code{prob}=\eqn{c(p_ {1},...,p_{r})}{c(p1,...,pr)} is the sum of \eqn{r} binomially distributed random variables with different \eqn{p_{i}}{pi} (and, in case, with different \eqn{n_{i}}{ni}):  

Z=\eqn{\sum_{i=1}^{r} Z_{i}}{\sum Zi}, \verb{} \eqn{Z} ~ \eqn{gbinom}(\code{size},\code{prob}), with \eqn{Z_{i}}{Zi} ~ \eqn{binom(n_{i},p_{i}),\ i=1,\dots ,r}{binom(ni,pi), i=1,...,r}.

Since the sum of Bernoulli distributed random variables is binomially distributed, \eqn{Z} can be also defined as: 

Z=\eqn{\sum_{i=1}^{r}\sum_{j=1}^{n_{i}}Z_{ij}}{\sum \sum Zij}, with \eqn{Z_{ij}}{Zij} ~ \eqn{binom(1,p_{i}),\ j=1,...,n_{i}}{binom(1,pi),  j=1,...,ni}.

The pmf is obtained by  an algorithm which is based on the convolution of Bernoulli distributions.
See the references below for further information.

The quantile is defined as the smallest value \eqn{x} such that \eqn{ F(x) \geq  p}{ F(x) \ge p} , where F is the cumulative distribution function. 

\code{rgbinom} uses the inversion method (see Devroye, 1986).

}
\value{
\code{dgbinom} gives the pmf, \code{pgbinom} gives the cdf, \code{qgbinom} gives the quantile function and \code{rgbinom} generates random deviates.
}
\references{
D.Kurz, H.Lewitschnig, J.Pilz, \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, IEEE Transactions on Reliability, Vol. 63, No. 2, June 2014.

K.J. Klauer, \emph{Kriteriumsorientierte Tests}, Verlag fuer Psychologie, Hogrefe, 1987, Goettingen, p. 208 ff.

M.Fisz, \emph{Wahrscheinlichkeitsrechnung und mathematische Statistik}, VEB Deutscher Verlag der Wissenschaften, 1973, p. 164 ff.

L.Devroye, Non-Uniform Random Variate Generation, Springer-Verlag, 1986, p. 85 ff.
}
\note{
If \code{size} contains just one trial number and \code{prob} one success probability, then the generalized binomial distribution results in the binomial distribution. 


The generalized binomial distribution described here is also known as Poisson-binomial distribution. See the link below to the package \code{poibin} for further information.
}

\seealso{
 \code{\link[poibin]{ppoibin}}, for another implementation of this distribution. \cr
\code{\link{dbinom}}
}

\examples{

## n=10 defect devices, divided in 3 failure types n1=2, n2=5, n3=3.
## 3 countermeasures with effectivities p1=0.8, p2=0.7, p3=0.3 are available.
## use dgbinom() to get the probabilities for x=0,...,10 failures solved.
dgbinom(x=c(0:10),size=c(2,5,3),prob=c(0.8,0.7,0.3))

## generation of N=100000 random values
rgbinom(100000,size=c(2,5,3),prob=c(0.8,0.7,0.3))


## n1=100, n2=100, n3=200, p1=0.001, p2=0.005, p3=0.01 
dgbinom(c(0:2),size=c(100,100,200),prob=c(0.001,0.005,0.01))
#  0.07343377 0.19260317 0.25173556
pgbinom(2,size=c(100,100,200),prob=c(0.001,0.005,0.01),lower.tail=FALSE)
# 0.4822275
}



