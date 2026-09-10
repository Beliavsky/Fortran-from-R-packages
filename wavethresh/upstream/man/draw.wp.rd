\name{draw.wp}
\alias{draw.wp}
\title{Draw wavelet packet associated with a wp object.}
\description{
This function draws a wavelet packet associated with a \code{\link{wp.object}}. 
}
\usage{
\method{draw}{wp}(wp, level, index, plot.it=TRUE, main, sub, xlab, ylab, \dots)

}
\arguments{
\item{wp}{The \code{\link{wp}} class object whose associated wavelet packet you wish to draw. }
\item{level}{The resolution level of wavelet packet in the wavelet packet decomposition that you wish to draw (corresponds to scale).}
\item{index}{The packet index of the wavelet packet in the wavelet packet decomposition that you wish to draw (corresponds to number of oscillations).}
\item{plot.it}{If TRUE then the wavelet packet is plotted on the active graphics device. If FALSE then the y-coordinates of the packet are returned. Note that x-coordinates are not returned (the packet is periodic on its range anyway). }
\item{main}{The main argument for the plot}
\item{sub}{The subtitle for the plot}
\item{xlab}{The labels for the x axis}
\item{ylab}{The labels for the y axis}
\item{\dots}{Additional arguments to pass to the \code{\link{drawwp.default}} function which does the drawing. In particular, arguments can be set to choose between drawing the mother wavelet and scaling function, to set the resolution of the plot, to choose between drawing one and two dimensional pictures. }
}
\details{
This function extracts the filter component from the \code{\link{wp}} object (which is constructed using the \code{\link{filter.select}} function) to decide which wavelet packet family to draw. Once decided the \code{\link{drawwp.default}} function is used to actually do the drawing. }
\value{
If the \code{plot.it} argument is set to \code{TRUE} then nothing is returned. Otherwise, if \code{plot.it} is set to \code{FALSE} the coordinates of what would have been plotted are returned. 
}
\note{If the \code{plot.it} argument is \code{TRUE} (which it is by default) a plot of the appropriate wavelet packet is plotted on the active graphics device.}
\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1998 }
\seealso{
\code{\link{filter.select}}, \code{\link{wp}}, \code{\link{wp.object}}, \code{\link{drawwp.default}}. 
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
#
# Now do the wavelet packet transform of the data using the Daubechies
# least-asymmetric wavelet N=10 (the default arguments in
# wp).
#
tdwp <- wp(test.data)
#
# What happens if we try to draw this new tdwp object?
#
\dontrun{draw(tdwd, level=4, index=12)}
}
\keyword{hplot}
\author{G P Nason}
