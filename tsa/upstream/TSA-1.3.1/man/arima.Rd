\name{arima}
\alias{arima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting an ARIMA model with Exogeneous Variables}
\description{
This function is identical to the arimax function which builds on and 
extends the capability of the arima function in R stats by allowing the
incorporation of transfer functions, and innovative and additive outliers.
For backward compatitibility, the function is also named arima.
Note in the computation of AIC, the number of parameters excludes the noise
variance. This function is heavily based on 
the arima function of the stats core
of R, see the help page of this function for details on arguments
\code{x} to \code{kappa}.}
\usage{
arima(x, order = c(0, 0, 0), seasonal = list(order = c(0, 0, 0), period = NA),
 xreg = NULL, include.mean = TRUE, transform.pars = TRUE, fixed = NULL, 
init = NULL, method = c("CSS-ML", "ML", "CSS"), n.cond, optim.control = list(),
 kappa = 1e+06, io = NULL, xtransf, transfer = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time series response}
  \item{order}{regular ARIMA order}
  \item{seasonal}{seasonal ARIMA order}
  \item{xreg}{a dataframe containing covariates}
  \item{include.mean}{if true, an intercept term is incorporated in the model;
applicable only to stationary models.}
  \item{transform.pars}{if true, the AR parameters are transformed to 
ensure stationarity}
  \item{fixed}{a vector indicating which coefficients are fixed or free}
  \item{init}{initial values}
  \item{method}{estimation method}
  \item{n.cond}{number of initial values to be conditioned on in a conditional
analysis}
  \item{optim.control}{control parameters for the optimization procedure}
  \item{kappa}{prior variance; used in dealing with initial values}
  \item{io}{a list of time points at which the model may have an innovative outlier. The time point of the outlier can be given either as absolute time point
 or as c(a,b), i.e. at the b-th 'month' of the a-th 'year' where each 
 year has frequency(x) months, assuming x is a time series.
}
  \item{xtransf}{
 xtranf is a matrix with each column containing a covariate that affects
 the time series response in terms of an ARMA filter of order (p,q), i.e.
 if Z is one such covariate, its effect on the time series is 
 \eqn{(theta_0+theta_1B+...+theta_{q-1}B^{q-1})/(1-phi_1 B -...-phi_p B^p) Z_t}
 In particular, if \eqn{p=0} and \eqn{q=1}, this specifies a simple regression 
relationship, which should be included in xreg and not here.
 Note that the filter starts with zero initial values. Hence, it is pertinent
 to mean-delete each distributed-lag covariate, and this 
is not done automatically.
}
  \item{transfer}{a list consisting of the ARMA orders for each transfer (distributed lag) covariate.  
}
}
\value{
An Arimax object contining the model fit.
}
\author{Original author of the arima function in R stats: Brian Ripley. The arimax function is based on the stats:::arima function, with modifications by Kung-Sik Chan.}
\seealso{ \code{\link{arima}}}
\examples{
data(hare)
arima(sqrt(hare),order=c(3,0,0))


}
\keyword{methods}
