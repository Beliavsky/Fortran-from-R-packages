\name{McLeod.Li.test}
\alias{McLeod.Li.test}
\title{McLeod-Li test} 
\description{
Perform the McLeod-Li test for conditional heteroscedascity
(ARCH).
}
\usage{
McLeod.Li.test(object, y, gof.lag, col = "red", omit.initial = TRUE, 
plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted Arima model, ususally the output from the arima
function. If supplied, then the Mcleod-Li test is applied to the residuals
of the model, and the y-argument is ignored.}
  \item{y}{time series data with which one wants to test for the presence of
conditional heteroscedascity}  
  \item{gof.lag}{maximum number of lags for which the test is carried out.}
  \item{col}{color of the reference line}
  \item{omit.initial}{suppress the initial (d+Ds) residuals if set to be TRUE}
  \item{plot}{suppress plotting if set to be FALSE}
  \item{\dots}{ other arguments to be passed to the plot function}
}
\details{
The test checks for the presence of conditional heteroscedascity by computing the
Ljung-Box (portmanteau) test with the squared data (if y is supplied and object
suppressed) or with the squared residuals from an arima model (if an arima
model is passed to the function via the object argument.) 
}
\value{
  \item{pvlaues}{the vector of p-values for the Ljung-Box test statistics
computed using the first \eqn{m} lags of the ACF of the squared data or residuals,
for \eqn{m} ranging from 1 to gof.lag.}
}
\references{
McLeod, A. I. and W. K. Li (1983). Diagnostic checking ARMA time series models
using squared residual autocorrelations. Journal of Time Series Analysis, 4,
269273.}
\author{Kung-Sik Chan}
\examples{
data(CREF)
r.cref=diff(log(CREF))*100
McLeod.Li.test(y=r.cref)
}
\keyword{methods}
