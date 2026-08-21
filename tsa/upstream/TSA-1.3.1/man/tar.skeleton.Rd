\name{tar.skeleton}
\alias{tar.skeleton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the asympotitc behavior of the skeleton of a TAR model}
\description{
The skeleton of a TAR model is obtained by suppressing the noise term from
the TAR model.  
}
\usage{
tar.skeleton(object, Phi1, Phi2, thd, d, p, ntransient = 500, n = 500, 
xstart, plot = TRUE,n.skeleton = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a TAR model fitted by the tar function; if it is supplied, the 
model parameters and initial values are extracted from it}
  \item{ntransient}{the burn-in size}
  \item{n}{sample size of the skeleton trajectory}
  \item{Phi1}{the coefficient vector of the lower-regime model}
  \item{Phi2}{the coefficient vector of the upper-regime model}
  \item{thd}{threshold}
  \item{d}{delay}
  \item{p}{maximum autoregressive order}
  \item{xstart}{initial values for the iteration of the skeleton}
  \item{plot}{if True, the time series plot of the skeleton is drawn} 
  \item{n.skeleton}{number of last n.skeleton points of the skeleton to be plotted}
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
\value{A vector that contains the trajectory of the skeleton, with the burn-in
discarded.}
\references{Tong, H. (1990) "Non-linear Time Series, a Dynamical System Approach," Clarendon Press Oxford.
"Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\seealso{ \code{\link{tar}}}
\examples{
data(prey.eq)
prey.tar.1=tar(y=log(prey.eq),p1=4,p2=4,d=3,a=.1,b=.9,print=TRUE)
tar.skeleton(prey.tar.1)
}
\keyword{methods}
