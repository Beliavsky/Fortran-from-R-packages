\name{draw.wst}
\alias{draw.wst}
\title{Draw mother wavelet or scaling function associated with wst object.}
\description{
This function draws the mother wavelet or scaling function associated with a  \code{\link{wst.object}}. 
}
\usage{
\method{draw}{wst}(wst, \dots)
}
\arguments{
\item{wst}{The \code{\link{wst}} class object whose associated wavelet or scaling function you wish to draw. }
\item{\dots}{Additional arguments to pass to the \code{\link{draw.default}} function which does the drawing. In particular, arguments can be set to choose between drawing the mother wavelet and scaling function, to set the resolution of the plot, to choose between drawing one and two dimensional pictures. }
}
\details{
This function extracts the \code{filter} component from the \code{\link{wst}} object (which is constructed using the \code{\link{filter.select}} function) to decide which wavelet packet family to draw. Once decided the \code{\link{draw.default}} function is used to actually do the drawing. }
\value{
If the \code{plot.it} argument is set to \code{TRUE} then nothing is returned. Otherwise, Otherwise, as with \code{\link{draw.default}}, the coordinates of what would have been plotted are returned. 
}
\note{If the \code{plot.it} argument is \code{TRUE} (which it is by default) a plot of the appropriate wavelet packet is plotted on the active graphics device.}
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1995 }
\seealso{
\code{\link{filter.select}}, \code{\link{wst.object}}, \code{\link{draw.default}}. 
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
#
# Now do the \code{packet-ordered non-decimated DWT}  of the data using the Daubechies
# least-asymmetric wavelet N=10 (the default arguments in \code{\link{wst}}).
#
tdwst <- wst(test.data)
#
# What happens if we try to draw this new tdwst object?
#
\dontrun{draw(tdwst)}
#
# We get a picture of the wavelet that did the transform
#
}
\keyword{hplot}
\author{G P Nason}
