\name{garch.sim}
\alias{garch.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate a GARCH process}
\description{Simulate a GARCH process.}
\usage{
garch.sim(alpha, beta, n = 100, rnd = rnorm, ntrans = 100,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{The vector of  ARCH coefficients including the intercept term as the 
first element}
  \item{beta}{The vector of GARCH coefficients}
  \item{n}{sample size}
  \item{rnd}{random number generator for the noise; default is normal}
  \item{ntrans}{burn-in size, i.e. number of initial simulated data to be discarded}
\item{...}{parameters to be passed to the random number generator}
}
\details{Simulate data from the GARCH(p,q) model:
\eqn{x_t=\sigma_{t|t-1} e_t} where \eqn{\{e_t\}} is iid,  \eqn{e_t} 
independent of past \eqn{x_{t-s}, s=1,2,\ldots}, and
\deqn{\sigma_{t|t-1}=\sum_{j=1}^p \beta_j \sigma_{t-j|t-j-1}+
\alpha_0+\sum_{j=1}^q \alpha_j x_{t-i}^2}
}
\value{simulated GARCH time series of size n.
}
\author{Kung-Sik Chan}
\examples{
set.seed(1235678)
garch01.sim=garch.sim(alpha=c(.01,.9),n=500)
plot(garch01.sim,type='l', main='',ylab=expression(r[t]),xlab='t')
}
\keyword{methods}
