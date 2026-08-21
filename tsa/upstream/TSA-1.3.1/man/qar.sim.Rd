\name{qar.sim}
\alias{qar.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate a first-order quadratic AR model}
\description{
Simulates a first-order quadratic AR model with normally distributed noise.
}
\usage{
qar.sim(const = 0, phi0 = 0, phi1 = 0.5, sigma = 1, n = 20, init = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{const}{intercept}
  \item{phi0}{coefficient of the lag 1}
  \item{phi1}{coefficient of the squared lag 1}
  \item{sigma}{noise standard deviation}
  \item{n}{sample size}
  \item{init}{number of burn-in values}
}
\details{
The quadratic AR(1) model specifies that 
\deqn{Y_t = \mathrm{const}+\phi_0 Y_{t-1}+\phi_1 Y^2_{t-1}+e_t}
where \eqn{e_t} are iid normally distributed with zero mean and standard
deviation \eqn{\sigma}. If \eqn{\sigma=0}, the model is deterministic. 
}
\value{A simulated series from the quadratic AR(1) model, as a vector}  
\author{Kung-Sik Chan}
\seealso{ \code{\link{tar.sim}}}
\examples{
set.seed(1234567)
plot(y=qar.sim(n=15,phi1=.5,sigma=1),x=1:15,type='l',ylab=expression(Y[t]),xlab='t')
y=qar.sim(n=100,const=0.0,phi0=3.97, phi1=-3.97,sigma=0,init=.377)
plot(y,x=1:100,type='l',ylab=expression(Y[t]),xlab='t')
acf(y,main='')
}
\keyword{methods}
