\name{plot.wst2D}
\alias{plot.wst2D}
\title{Plot packet-ordered 2D non-decimated wavelet coefficients.}
\description{
This function plots packet-ordered 2D non-decimated wavelet coefficients arising from a \code{\link{wst2D}} object. 
}
\usage{
\method{plot}{wst2D}(x, plot.type="level", main="", ...)
}
\arguments{
\item{x}{The \code{\link{wst2D}} object whose coefficients you wish to plot.}
\item{plot.type}{So far the only valid argument is "level" which plots coefficients a level at a time.}
\item{main}{The main title of the plot.} 
\item{...}{Any other arguments.}

}
\details{
The coefficients in a \code{\link{wst2D}} object are stored in a three-dimensional subarray called \code{wst2D}. The first index of the 3D array indexes the resolution level of coefficients: this function with \code{plot.type="level"} causes an image of coefficients to be plotted one for each resolution level. 

The following corresponds to images produced on S+ graphics devices (e.g. image on \code{motif()}). Given a resolution level there are \code{4^(nlevelsWT-level)} packets within a level. Each packet can be addressed by a base-4 string of length \code{nlevels-level}. A zero corresponds to no shift, a 1 to a horizontal shift, a 2 to a vertical shift and a 3 to both a horizontal and vertical shift. 

So, for examples, at resolution level \code{nlevelsWT-1} there are 4 sub-images each containing 4 sub-images. The main subimages correspond to (clockwise from bottom-left) no shift, horizontal shift, both shift and vertical shifts. The sub-images of the sub-images correspond to the usual smooth, horizontal detail, diagonal detail and vertical detail (clockwise, again from bottom left). Coarser resolution levels correspond to finer shifts! The following figure demonstrates the \code{nlevels-1} resolution level for the \code{ua} image (although the whole image has been rotated by 90 degrees clockwise for display here!):
}
\value{
A plot of the coefficients contained within the \code{\link{wst2D}} object is produced. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{getpacket.wst2D}}, \code{\link{putpacket.wst2D}}, \code{\link{wst2D}}, \code{\link{wst2D.object}}.  
}
\examples{
#
# The above picture is one of a series produced by 
#
#plot(uawst2D)
#
# Where the uawst2D object was produced in the EXAMPLES section
# of the help for \code{\link{wst2D}}
}
\keyword{hplot}
\author{G P Nason}
