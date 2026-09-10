\name{draw.imwd}
\alias{draw.imwd}
\title{Draw mother wavelet associated with an imwd object. }
\description{
This function draws the mother wavelet associated with an \code{\link{imwd.object}} --- a two-dimensional wavelet decomposition. 
}
\usage{
\method{draw}{imwd}(wd, resolution=128, ...)
}
\arguments{
\item{wd}{The \code{\link{imwd}} class object whose associated wavelet you wish to draw. }
\item{resolution}{The resolution at which the computation is done to compute the wavelet picture. Generally the resolution should be lower for two-dimensional wavelets since the number of computations is proportional to the square of the resolution (the DWT is still O(n) though).}
\item{\dots}{Additional arguments to pass to the \code{\link{draw.default}} function which does the drawing.}
}
\details{
This function extracts the \code{filter} component from the \code{\link{imwd}} object (which is constructed using the \code{\link{filter.select}} function) to decide which wavelet to draw. Once decided the \code{\link{draw.default}} function is used to actually do the drawing. 
}
\value{
If the plot.it argument is set to \code{TRUE} then nothing is returned. Otherwise, as with \code{\link{draw.default}}, the coordinates of what would have been plotted are returned. 
}
\note{
If the \code{plot.it} argument is\code{TRUE} (which it is by default) a plot of the mother wavelet or scaling function is plotted on the active graphics device.
}
\section{RELEASE}{Version 2 Copyright Guy Nason 1993 }
\seealso{
\code{\link{filter.select}}, \code{\link{imwd.object}}, \code{\link{draw.default}}. 
}
\examples{
#
# Let's use the lennon test image
#
data(lennon)
\dontrun{image(lennon)}
#
# Now let's do the 2D discrete wavelet transform using Daubechies'
# least-asymmetric wavelet N=6
#
lwd <- imwd(lennon, filter.number=6)
#
# And now draw the wavelet that did this transform
#
\dontrun{draw(lwd)}
#
# A nice little two-dimensional wavelet!
#

}
\keyword{hplot}
\author{G P Nason}
