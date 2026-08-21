\name{tsdiag.TAR}
\alias{tsdiag.TAR}
\title{Model diagnostics for a fitted TAR model}
\description{
The time series plot and the sample ACF of the standardized residuals are plotted.
Also, a portmanteau test for detecting residual correlations in the standardized
residuals are carried out. 
}
\usage{
\method{tsdiag}{TAR}(object, gof.lag, col = "red",xlab = "t", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted TAR model output from the tar function}
  \item{gof.lag}{number of lags of ACF to be examined}
  \item{col}{color of the lines flagging outliers, etc.}
  \item{xlab}{x labels for the plots}
  \item{\dots}{ any additional user-supplied options to be passed to the
acf function}
}
\value{
Side effects: plot the time-series plot of the standardized residuals, their
sample ACF and portmanteau test for residual autocorrelations in the standardized errors.} 
\references{"Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\seealso{ \code{\link{tar}}}
\examples{
data(prey.eq)
prey.tar.1=tar(y=log(prey.eq),p1=4,p2=4,d=3,a=.1,b=.9,print=TRUE)
tsdiag(prey.tar.1)
}
\keyword{methods}
