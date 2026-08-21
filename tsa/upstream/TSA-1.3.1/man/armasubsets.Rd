\name{armasubsets}
\alias{armasubsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of Subset ARMA Models} 
\description{
This function finds a number of subset ARMA models.
 A "long" AR model is fitted to the data y to compute the residuals
 which are taken as a proxy of the error process.
 Then, an ARMA model is approximated by a regression model with the 
 the covariates being the lags of the time series and the lags of the
 error process.
 Subset ARMA models may then be selected using the subset regression 
 technique by leaps and bounds, via the regsubsets function of the leaps
 package in R.
}
\usage{
armasubsets(y, nar, nma, y.name = "Y", ar.method = "ols", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ \item{y}{time-series data}
  \item{nar}{maximum AR order}
  \item{nma}{maximum MA order}
  \item{y.name}{label of the time series}
  \item{ar.method}{method used for fitting the long AR model; default is ols
with the AR order determined by AIC}
  \item{\dots}{arguments passed to the plot.armasubsets function} }
\value{
An object of the armasubsets class to be processed by the plot.armasubsets
function.} 
\author{Kung-Sik Chan}
\examples{
set.seed(92397)
test=arima.sim(model=list(ar=c(rep(0,11),.8),ma=c(rep(0,11),0.7)),n=120)
res=armasubsets(y=test,nar=14,nma=14,y.name='test',ar.method='ols')
plot(res)
}
\keyword{methods}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
