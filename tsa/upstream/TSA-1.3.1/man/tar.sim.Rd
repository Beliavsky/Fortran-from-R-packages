\name{tar.sim}
\alias{tar.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate a two-regime TAR model}
\description{Simulate a two-regime TAR model.}
\usage{
tar.sim(object, ntransient = 500, n = 500, Phi1, Phi2, thd, d, p, sigma1, 
sigma2, xstart = rep(0, max(p,d)), e)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a TAR model fitted by the tar function; if it is supplied, the 
model parameters and initial values are extracted from it}
  \item{ntransient}{the burn-in size}
  \item{n}{sample size of the simulated series}
  \item{Phi1}{the coefficient vector of the lower-regime model}
  \item{Phi2}{the coefficient vector of the upper-regime model}
  \item{thd}{threshold}
  \item{d}{delay}
  \item{p}{maximum autoregressive order}
  \item{sigma1}{noise std. dev. in the lower regime}
  \item{sigma2}{noise std. dev. in the upper regime}
  \item{xstart}{initial values for the simulation}
  \item{e}{standardized noise series of size equal to length(xstart)+ntransient+n; if missing, it will be generated as some normally distributed errors}
}
\details{
The two-regime Threshold Autoregressive  (TAR) model is given by the following
formula:
\deqn{
        Y_t = \phi_{1,0}+\phi_{1,1} Y_{t-1} +\ldots+                                                                                          \phi_{1,p} Y_{t-p_1} +\sigma_1 e_t,
         \mbox{  if } Y_{t-d}\le r }
\deqn{ Y_t = \phi_{2,0}+\phi_{2,1} Y_{t-1} +\ldots+\phi_{2,p_2} Y_{t-p}+\sigma_2 e_t,
\mbox{  if } Y_{t-d} > r.}
where r is the threshold and d the delay.
}
\value{A list containing the following components:
  \item{y}{simulated TAR series}
  \item{e}{the standardized errors} 
  ...
}
\references{Tong, H. (1990) "Non-linear Time Series, a Dynamical System Approach," Clarendon Press Oxford
"Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\seealso{ \code{\link{tar}}}
\examples{
set.seed(1234579)
y=tar.sim(n=100,Phi1=c(0,0.5),
Phi2=c(0,-1.8),p=1,d=1,sigma1=1,thd=-1,
sigma2=2)$y
plot(y=y,x=1:100,type='b',xlab="t",ylab=expression(Y[t]))
}
\keyword{methods}
