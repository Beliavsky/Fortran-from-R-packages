\name{tacvfARMA}
\alias{tacvfARMA}
\title{ theoretical autocovariance function (acvf) of ARMA }
\description{
The theoretical autocovariance function of ARMA(p,q) process is computed.
This is more useful in some situations than the built-in R function \code{\link{ARMAacf}}.
See Details.
}

\usage{
tacvfARMA(phi = numeric(0), theta = numeric(0), maxLag = 1, sigma2 = 1)
}

\arguments{
  \item{phi}{ ar parameters }
  \item{theta}{ ma parameters }
  \item{maxLag}{ acvf is computed at lags 0, ..., maxLag }
  \item{sigma2}{ innovation variance }
}

\details{
The details of the autocovariance computation are given in McLeod (1975).

In addition to this computation, we also test if the model is stationary-causal
or not.  The test, which is included directly in the function, uses
the Durbin-Levison recursion to transform from the phi parameters to the pacf.
See McLeod and Zhang (2006, eqn. (1)) for more details.
Formally, the stationary-causal condition requires that all roots of the
polynomial equation, 
\deqn{1 - phi[1]*B -...- phi[p]*B^p = 0 } 
must lie outside
the unit circle (Brockwell and Davis, 1991, Section 3.3).

This function is included because it is necessary to demonstrate that in
the case of ARMA models, TrenchInverse and the built-in R function
predict.Arima produce equivalent results.
See Example 1 in the documentation for \code{\link{TrenchForecast}}
and the example discussed in McLeod, Yu and Krougly (2007, 3.2). 
}

\value{
Vector of length maxLag containing the autocovariances at
lags 0, ..., maxLag. But see Warning below.
}

\references{  
P.J. Brockwell and R.A. Davis (1991)
Time Series: Theory and Methods. Springer.

A.I. McLeod (1975) 
Derivation of the theoretical autocovariance function of autoregressive-moving
average models, Applied Statistics 24, 255-256.

A.I. McLeod and Zhang, Y. (2006)
Partial autocorrelation parameterizations for subset autoregression,
Journal of Time Series Analysis,  

McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.
}

\note{
An error is returned if the model is not stationary-causal.
}

\author{A.I. McLeod }

\seealso{ 
\code{\link{ARMAacf}}
}

\examples{
#Example 1.  Estimate the acvf of a fitted ARMA model
#There are two methods but they give slighly different results,
#general script, just change z, p, q, ML
z<-sqrt(sunspot.year)
n<-length(z)
p<-9
q<-0
ML<-5
#for different data/model just reset above
out<-arima(z, order=c(p,0,q))
phi<-theta<-numeric(0)
if (p>0) phi<-coef(out)[1:p]
if (q>0) theta<-coef(out)[(p+1):(p+q)]
zm<-coef(out)[p+q+1]
sigma2<-out$sigma2
rA<-tacvfARMA(phi, theta, maxLag=n+ML-1, sigma2=sigma2)
rB<-var(z)*ARMAacf(ar=phi, ma=theta, lag.max=n+ML-1)
#rA and rB are slighly different
cbind(rA[1:5],rB[1:5])
#
#Example 2. Compute Rsq for fitted ARMA model
#Rsq = 1 - (series variance / innovation variance)
#Again there are two methods but only the first method is guaranteed to
#produce an Rsq which is non-negative!
#Run last example and then evaluate the script below:
RsqA <- 1 - rA/sigma2
RsqB <- 1 - rB/sigma2
#
#Example 3. Test if model is stationary-causal or not.
StationaryQ <- function(phi) tryCatch(is.vector(tacvfARMA(phi=phi)),error=function(e) FALSE )
StationaryQ(1.1) #AR(1) with phi=1.1 is not stationary-causal.
#try with parameters from Example 1 above
StationaryQ(phi)
}

\keyword{ ts }
