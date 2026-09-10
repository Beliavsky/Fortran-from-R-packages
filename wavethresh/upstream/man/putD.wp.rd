\name{putD.wp}
\alias{putD.wp}
\title{Puts a whole resolution level of wavelet packet coeffients into wp wavelet object.}
\description{
Makes a copy of the \code{\link{wp}} object, replaces a whole resolution level of wavelet packet coefficients data in the copy, and then returns the copy. 
}
\usage{
\method{putD}{wp}(wp, level, value, \dots)
}
\arguments{
\item{wp}{Wavelet packet object into which you wish to insert the wavelet packet coefficients.}
\item{level}{the resolution level at which you wish to replace the wavelet packet coefficients.} 
\item{value}{the replacement data, this should be of the correct length.} 
\item{\dots}{any other arguments}
}
\details{
The function \code{\link{accessD.wp}} obtains the wavelet packet coefficients for a particular level.
 
For wavelet packet transforms the number of coefficients at each resolution level is the same and equal to \code{2^nlevelsWT} where \code{nlevels} is the number of levels in the \code{\link{wp.object}}. The number of coefficients at each resolution level is also, of course, the number of data points used to initially form the \code{\link{wp}} object in the first place. 

Use the \code{\link{accessD.wp}} to extract whole resolution levels of wavelet packet coefficients. 

We don't recommend that you use this function unless you really know what you are doing. Usually it is more convenient to manipulate individual \emph{packets} of coefficients using \code{\link{getpacket}}/\code{\link{putpacket}} functions. If you must use this function to insert whole resolution levels of coefficients you must ensure that the data vector you supply is valid: i.e. contains packet coefficients in the right order. 
}
\value{
A \code{\link{wp}} class object containing the modified wavelet packet coefficients. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wp.object}}, \code{\link{wp}}, \code{\link{accessD}}, \code{\link{accessD.wp}}, \code{\link{getpacket.wp}}, \code{\link{putpacket.wp}}. 
}
\examples{
#
# Generate an EMPTY wp object:
#
zero <- rep(0, 16)
zerowp <- wp(zero)
#
# Put some random mother wavelet coefficients into the object at
# resolution level 2. For the wavelet packet transform there
# are always 16 coefficients at every resolution level. 
#
mod.zerowp <- putD( zerowp, level=2, v=rnorm(16))
#
# If you plot mod.zerowp you will see that there are only 
# coefficients at resolution level 2 where you just put the coefficients.
}
\keyword{manip}
\author{G P Nason}
