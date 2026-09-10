\name{denplot}
\alias{denplot}
\title{Calculate plotting information for a density estimate.}
\usage{
denplot(wr, coef, nT=20, lims, n=50)
}
\arguments{
\item{wr}{Scaling function coefficients, usually at some high level and
	usually smoothed (thresholded).}
\item{coef}{The output from \code{\link{denproj}} for this analysis,
i.e. the object containing the empirical scaling function coefficients.
This is required because of the information it contains about the wavelet
filter used, the resolution of the projection,
and the bounds on the translation index of the scaling function coefficients.}
\item{lims}{Vector containing the minimum and maximum x values required on the
	plot.}
\item{nT}{The number of iterations to be performed in the
	Daubechies-Lagarias algorithm,
	which is used to evaluate the scaling functions of the specified
	wavelet basis at the plotting points.}
\item{n}{The number of points at which the density estimate is to be evaluated.}
}
\description{
Calculates plotting information for a wavelet density estimate from high level
scaling function coefficients.
}
\details{
The density estimate is evaluated at \code{n} points between the values
in \code{lims}. This function can be used to plot the empirical
scaling function density estimate by entering
\code{wr=coef$coef}, but since the empirical coefficients are usually found
at some very high resolution,
such a plot will be very noisy and not very informative.
This function will be of much more use as and when
thresholding function are included in this density estimation package.
}
\value{
A list with components:
\item{x}{The points at which the density estimate is evaluated.}
\item{y}{The values of the density estimate at the points in \code{x}.}
}
\examples{
# Simulate data from the claw density and find the 
# empirical scaling function coefficients at a lowish resolution and plot
# the resulting density estimate
data <- rclaw(100)
datahr <- denproj(data, J=3, filter.number=2,family="DaubExPhase")
datapl <- denplot(datahr$coef, datahr, lims=c(-3,3), n=1000)
\dontrun{plot(datapl, type="l")}
}
\seealso{\code{\link{denproj}},\code{\link{rclaw}}}
\author{David Herrick}
\keyword{smooth}
