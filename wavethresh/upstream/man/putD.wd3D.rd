\name{putD.wd3D}
\alias{putD.wd3D}
\title{Put wavelet coefficient array into a 3D wavelet object}
\description{
This function put an array of wavelet coefficients, corresponding to a particular resolution level into a \code{\link{wd}} wavelet decomposition object. 

The pyramid of coefficients in a wavelet decomposition (returned from the \code{\link{wd3D}} function, say) are packed into a single array in \code{WaveThresh3}.
}
\usage{
\method{putD}{wd3D}(x, v, \dots)
}
\arguments{
\item{x}{3D Wavelet decomposition object into which you wish to insert the wavelet coefficients.}
\item{v}{This argument is a list with the following components: 
	\describe{
     \item{a}{A 3-dimensional array with each dimension of length equal to two to the power of lev which is the level at which you wish to insert the coefficients into x.}
     \item{lev}{The level at which you wish to insert the coefficients into \code{x}.} 
     \item{block}{A character string indicating which coefficient block you wish to insert the coefficients into. This can be one of GGG, GGH, GHG, GHH, HGG, HGH, HHG. Additionally this can be HHH when the lev argument above is zero.}}}
\item{\dots}{Other arguments}
}
\details{
The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as an array. 

Note that this function is a method for the generic function \code{\link{putD}}. 
}
\value{
A new \code{\link{wd3D.object}} is returned with the coefficients at level \code{lev} in block given by block are replaced by the contents of \code{a}, if \code{a} is of the correct dimensions! 
}
\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{accessD}}, \code{\link{accessD.wd3D}}, \code{\link{print.wd3D}}, \code{\link{putD}}, \code{\link{putDwd3Dcheck}}, \code{\link{summary.wd3D}}, \code{\link{threshold.wd3D}}, \code{\link{wd3D}}, \code{\link{wd3D.object}}, \code{\link{wr3D}}. 
}
\examples{
#
# Generate some test data
#
a <- array(rnorm(8*8*8), dim=c(8,8,8))
#
# Perform the 3D DWT
#
awd3D <- wd3D(a)
#
# Replace the second level coefficients by uniform random variables
# in block GGG (for some reason)
#
#
newsubarray <- list(a = array(runif(4*4*4), dim=c(4,4,4)), lev=2, block="GGG")
awd3D <- putD(awd3D,  v=newsubarray)
}
\keyword{manip}
\author{G P Nason}
