\name{CWavDE}
\alias{CWavDE}
\title{Simple wavelet density estimator with hard thresholding}
\usage{
CWavDE(x, Jmax, threshold=0, nout=100, primary.resolution=1, filter.number=10,
	family="DaubLeAsymm", verbose=0, SF=NULL, WV=NULL)
}
\arguments{
\item{x}{Vector of real numbers. This is the data for which you want
	a density estimate for}
\item{Jmax}{The maximum resolution of wavelets}
\item{threshold}{The hard threshold value for the wavelet coefficients}
\item{nout}{The number of ordinates in the density estimate}
\item{primary.resolution}{The usual wavelet density estimator primary
	resolution}
\item{filter.number}{The wavelet filter number, see \code{\link{filter.select}}}
\item{family}{The wavelet family, see \code{\link{filter.select}}}
\item{verbose}{The level of reporting performed by the function, legit
	values are 0, 1 or 2, with 2 being more reports}
\item{SF}{Scaling function values in format as returned by
	\code{\link{draw.default}}}
\item{WV}{Wavelet function values in format as returned by
	\code{\link{draw.default}}}
}
\description{
This function implements the density estimator with hard thresholding
described by Hall, P. and Patil, P. (1995) Formulae for mean integrated squared
error of nonlinear wavelet-based density estimators,
\emph{Ann. Statist.}, \bold{23}, 905-928.
}
\details{
As the description.
}
\value{
A list containing the following components:
\item{x}{A vector of length \code{nout} that covers the range of the
	input data \code{x}, plus some more depending on the
	support of the wavelet and the primary resolution.}
\item{y}{A vector of length \code{nout} that contains the output
	wavelet density estimate}
\item{sfix}{The integer values of the translates of the scaling functions
	used in the estimate}
\item{wvixmin}{As for sfix, but a vector of length \code{Jmax} which
	contains the minimum integer wavelet translates}
\item{wvixmax}{As for wvixmin, but with the maxima}
}
\examples{
#
# Let's generate a bi-modal artificial set of data.
#
x <- c( rnorm(100), rnorm(100, 10))
#
# Now perform simple wavelet density estimate
#
wde <- CWavDE(x, Jmax=10, threshold=1)
#
# Plot results
#
\dontrun{plot(wde$x, wde$y, type="l")}
}
\author{G P Nason}
\keyword{smooth}
