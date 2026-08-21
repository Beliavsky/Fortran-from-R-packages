\name{periodogram}
\alias{periodogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computing the periodogram}
\description{
This is a wrapper that computes the periodogram}
\usage{
periodogram(y,log='no',plot=TRUE,ylab="Periodogram", xlab="Frequency",lwd=2,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{ \item{y}{A univariate time series}
  \item{log}{if set to "yes", the periodogram is plotted on the log-scale; default="no"}
  \item{plot}{The periodogram is plotted if it is set to be TRUE which is the
default}
  \item{ylab}{label on the y-axis}
  \item{xlab}{label on the x-axis}
  \item{lwd}{thickness of the periodogram lines}
  \item{...}{other arguments to be passed to the plot function} }
\value{
A list that contains the following elements:
\item{freq}{Vector of frequencies at which the spectral density is estimated. (Possibly approximate Fourier frequencies.} 
\item{spec}{Vector of estimates of the periodogram at frequencies corresponding to freq.} 
}
\references{Bloomfield, P. (1976) Fourier Analysis of Time Series: An Introduction. Wiley. 

Brockwell, P. J. and Davis, R. A. (1991) Time Series: Theory and Methods. Second edition. Springer. 
}
\examples{
data(star)
plot(star,xlab='Day',ylab='Brightness')
periodogram(star,ylab='Variable Star Periodogram');  abline(h=0)
}
\keyword{methods}
