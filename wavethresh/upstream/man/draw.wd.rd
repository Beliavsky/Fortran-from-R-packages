\name{draw.wd}
\alias{draw.wd}
\title{Draw mother wavelet or scaling function associated with wd object. }
\description{
This function draws the mother wavelet or scaling function associated with a \code{\link{wd.object}}. 
}
\usage{
\method{draw}{wd}(wd, ...)
}
\arguments{
\item{wd}{The \code{\link{wd}} class object whose associated wavelet or scaling function you wish to draw. }
\item{\dots}{Additional arguments to pass to the \code{\link{draw.default}} function which does the drawing. In particular, arguments can be set to choose between drawing the mother wavelet and scaling function, to set the resolution of the plot, to choose between drawing one and two dimensional pictures.}
}
\details{
This function extracts the filter component from the \code{\link{wd}} object (which is constructed using the \code{\link{filter.select}} function) to decide which wavelet to draw. Once decided the \code{\link{draw.default}} function is used to actually do the drawing. 
}
\value{
If the \code{plot.it} argument is set to TRUE then nothing is returned. Otherwise, as with \code{\link{draw.default}}, the coordinates of what would have been plotted are returned. 
}
\note{If the \code{plot.it} argument is \code{TRUE} (which it is by default) a plot of the mother wavelet or scaling function is plotted on the active graphics device.}
\section{RELEASE}{Version 2 Copyright Guy Nason 1993 }
\seealso{
\code{\link{filter.select}}, \code{\link{wd.object}}, \code{\link{draw.default}}. 
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
#
# Now do the discrete wavelet transform of the data using the Daubechies
# least-asymmetric wavelet N=10 (the default arguments in
# wd).
#
tdwd <- wd(test.data)
#
# What happens if we try to draw this new tdwd object?
#
\dontrun{draw(tdwd)}
#
# We get a picture of the wavelet that did the transform
#

}
\keyword{hplot}
\author{G P Nason}
