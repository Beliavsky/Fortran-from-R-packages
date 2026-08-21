\name{tar}
\alias{tar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of a TAR model}
\description{
Estimation of a two-regime TAR model. 
}
\usage{
tar(y, p1, p2, d, is.constant1 = TRUE, is.constant2 = TRUE, transform = "no",
 center = FALSE, standard = FALSE, estimate.thd = TRUE, threshold, 
method = c("MAIC", "CLS")[1], a = 0.05, b = 0.95, order.select = TRUE, print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{time series}
  \item{p1}{AR order of the lower regime}
  \item{p2}{AR order of the upper regime}
  \item{d}{delay parameter}
  \item{is.constant1}{if True, intercept included in the lower regime, otherwise
the intercept is fixed at zero} 
  \item{is.constant2}{similar to is.constant1 but for the upper regime}
  \item{transform}{available transformations: "no" (i.e. use raw data), "log", "log10" and
"sqrt"}
  \item{center}{if set to be True, data are centered before analysis}
  \item{standard}{if set to be True, data are standardized before analysis}
  \item{estimate.thd}{if True, threshold parameter is estimated, otherwise
it is fixed at the value supplied by threshold}
  \item{threshold}{known threshold value, only needed to be supplied if estimate.thd is set to be False.}
  \item{method}{"MAIC": estimate the TAR model by minimizing  the AIC; 
"CLS": estimate the TAR model by the method of Conditional Least Squares.}
  \item{a}{lower percent; the threshold is searched over the interval defined by the
a*100 percentile to the b*100 percentile of the time-series variable}
  \item{b}{upper percent}
  \item{order.select}{If method is "MAIC", setting order.select to True will
enable the function to further select the AR order in 
each regime by minimizing 
AIC}
  \item{print}{if True, the estimated model will be printed}
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
\value{
A list of class "TAR" which can be further processed by the 
by the predict and tsdiag functions.
}
\references{Tong, H. (1990) "Non-linear Time Series, a Dynamical System Approach," Clarendon Press Oxford
                                                                                
"Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\seealso{\code{\link{predict.TAR}}, 
\code{\link{tsdiag.TAR}}, 
\code{\link{tar.sim}}, 
\code{\link{tar.skeleton}} 
}
\examples{
data(prey.eq)
prey.tar.1=tar(y=log(prey.eq),p1=4,p2=4,d=3,a=.1,b=.9,print=TRUE)
}
\keyword{methods}
