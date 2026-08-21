\name{acf}
\alias{acf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auto- and Cross- Covariance and -Correlation Function Estimation}
\description{
This function calls the acf function in the stats package and processes to drop lag-0 of the acf. It only works for univariate time series, so x below should be 1-dimensional.
}
\usage{
acf(x, lag.max = NULL, type = c("correlation", "covariance", "partial")[1], 
plot = TRUE, na.action = na.fail, demean = TRUE, drop.lag.0 = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a univariate or multivariate (not ccf) numeric time series object or a numeric vector or matrix, or an "acf" object.}
  \item{lag.max}{maximum number of lags at which to calculate the acf. Default is 10*log10(N/m) where N is the number of observations and m the number of series.}
  \item{type}{character string giving the type of acf to be computed. Allowed values are "correlation" (the default), "covariance" or "partial".}
  \item{plot}{
logical. If TRUE (the default) the acf is plotted.}                                                                           
  \item{na.action}{function to be called to handle missing values. na.pass can be used.}
  \item{demean}{
logical. Should the covariances be about the sample means?}
  \item{drop.lag.0}{logical. Should lag 0 be dropped}
  \item{\dots}{further arguments to be passed to plot.acf.}
}
\value{
An object of class "acf", which is a list with the following elements: 

\item{lag}{ A three dimensional array containing the lags at which the acf is estimated.} 
\item{acf}{ An array with the same dimensions as lag containing the estimated acf.} 
\item{type}{ The type of correlation (same as the type argument).} 
\item{n.used}{ The number of observations in the time series.} 
\item{series}{ The name of the series x.} 
\item{snames}{ The series names for a multivariate time series.} 
}
\references{ ~put references to the literature/web site here ~ }
\author{Original authors of stats:::acf are: Paul Gilbert, Martyn Plummer, B.D. Ripley. This wrapper is written by Kung-Sik Chan} 
\seealso{\code{\link{plot.acf}}, \code{\link{ARMAacf}} for the exact autocorrelations of a given ARMA process.} 
\examples{
data(rwalk)
model1=lm(rwalk~time(rwalk))
summary(model1)
acf(rstudent(model1),main='')
}
\keyword{methods}
