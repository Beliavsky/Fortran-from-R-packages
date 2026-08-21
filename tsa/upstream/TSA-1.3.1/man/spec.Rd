\name{spec}
\alias{spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computing the spectrum}
\description{
This is a wrapper that allows the user to invoke either the spec.pgram
function or the spec.ar function in the stats pacakge.  Note that the seasonal attribute of
the data, if it exists, will be removed, for our preferred way of presenting
the output.} 
\usage{ spec(x, taper = 0, detrend = FALSE, demean = TRUE, method = c("pgram", 
    "ar")[1], ci.plot = FALSE, ylim = range(c(lower.conf.band, upper.conf.band)), 
    ...) }
%- maybe also 'usage' for other objects documented here.
\arguments{A list that contains the following:
  \item{x}{A univariate or multivariate time series}
  \item{taper}{amount of taper; 0 is the default}
  \item{detrend}{logical; if True, the data are detrended; default is False}
  \item{demean}{logical; if True, the data are centered; default is True}
  \item{method}{String specifying the method used to estimate the spectral density. Allowed methods are "pgram" (the default) and "ar".}
  \item{ci.plot}{logical; if True, the 95\% confidence band will be plotted.}
  \item{ylim}{Plotting parameter vector specifying the minimum and maximum of the 
y-axis.}
  \item{\dots}{other arguments}
}
\value{
The output is from the spec.pgram function or spec.ar function, and 
the following description of the output is taken from the help manual
of the spec function in the stats package. 
An object of class "spec", which is a list containing at least the following components: 

\item{freq}{Vector of frequencies at which the spectral density is estimated. (Possibly approximate Fourier frequencies.) The units are the reciprocal of cycles per unit time (and not per observation spacing): see Details below.} 
\item{spec}{Vector (for univariate series) or matrix (for multivariate series) of estimates of the spectral density at frequencies corresponding to freq. 
coh NULL for univariate series. For multivariate time series, a matrix containing the squared coherency between different series. Column i + (j - 1) * (j - 2)/2 of coh contains the squared coherency between columns i and j of x, where i < j.} 
\item{phase}{NULL for univariate series. For multivariate time series a matrix containing the cross-spectrum phase between different series. The format is the same as coh.} 
\item{series}{The name of the time series.} 
\item{snames}{For multivariate input, the names of the component series.} 
\item{method}{The method used to calculate the spectrum.} 

The result is returned invisibly if plot is true.
}
\references{Bloomfield, P. (1976) Fourier Analysis of Time Series: An Introduction. Wiley. 

Brockwell, P. J. and Davis, R. A. (1991) Time Series: Theory and Methods. Second edition. Springer. 

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S-PLUS. Fourth edition. Springer. (Especially pages 3927.) 
}
\examples{
set.seed(271435); n=200; phi=-0.6
y=arima.sim(model=list(ar=phi),n=n)
k=kernel('daniell',m=15)
sp=spec(y,kernel=k,main='',sub='',xlab='Frequency',
ylab='Log(Smoothed Sample Spectrum)',ci.plot=TRUE,ci.col='black')
lines(sp$freq,ARMAspec(model=list(ar=phi),sp$freq,plot=FALSE)$spec,lty=4)
abline(h=0)
}
\keyword{methods}
