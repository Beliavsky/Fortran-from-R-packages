\name{gBox}
\alias{gBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized Portmanteau Tests for GARCH Models}
\description{
Perform a goodness-of-fit test for the GARCH model by checking whether the 
standardized residuals are iid based on the ACF of the absolute residuals
or squared residuals.}
\usage{gBox(model, lags = 1:20, x, method = c("squared", "absolute")[1], plot = TRUE)} 
\arguments{
  \item{model}{fitted model from the garch function of the tseries library}
  \item{lags}{a vector of maximum ACF lags to be  used in the test}
  \item{x}{time series data to which the GARCH model is fitted}
  \item{method}{"squared": test is based on squared residuals; "absolute": test is
based on absolute residuals}
  \item{plot}{logical variable, if TRUE, the p-values of the tests are plotted}
}
\value{
  \item{lags}{lags in the input}
  \item{pvalue}{a vector of p-values of the tests}
  \item{method}{method used}
  \item{x}{x}
}
\references{"Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\examples{
\donttest{
require(tseries) # need to uncomment this line when running the example
data(CREF)
r.cref=diff(log(CREF))*100
m1=tseries::garch(x=r.cref,order=c(1,1))
summary(m1)
gBox(m1,x=r.cref,method='squared')
}
}
\keyword{methods}
