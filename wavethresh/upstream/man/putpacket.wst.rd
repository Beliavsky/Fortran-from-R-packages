\name{putpacket.wst}
\alias{putpacket.wst}
\title{Put a packet of coefficients into a packet ordered non-decimated wavelet object (wst).}
\description{
This function inserts a packet of coefficients into a packet-ordered non-decimated wavelet object (\code{\link{wst}}) object. The \code{\link{wst}} objects are computed by the \code{\link{wst}} function amongst others. 
}
\usage{
\method{putpacket}{wst}(wst, level, index, packet, \dots)
}
\arguments{
\item{wst}{Packet-ordered non-decimated wavelet object into which you wish to insert the packet.} 
\item{level}{The resolution level of the coefficients that you wish to insert.} 
\item{index}{The index number within the resolution level of the packet of coefficients that you wish to insert.}
\item{packet}{A vector of coefficients that you wish to insert into the \code{\link{wst}} object. The length that the packet has to may be determined by extracting the same packet of coefficients using the \code{\link{getpacket.wst}} function and using the S-Plus length function to determine the length!} 
\item{\dots}{any other arguments}
}
\details{
This function actually calls the \code{\link{putpacket.wp}} to do the insertion.
 
In the future this function will be extended to insert father wavelet coefficients as well. 
}
\value{
An object of class \code{\link{wst.object}} containing the packet ordered non-decimated wavelet coefficients that have been modified: i.e. with packet inserted.
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{getpacket.wst}}, \code{\link{putpacket}}, \code{\link{putpacket.wp}}, \code{\link{wst}}, \code{\link{wst.object}}.
 }
\examples{
#
# Take the packet-ordered non-decimated transform of some random data 
#
MyWST <- wst(rnorm(1:512))
#
# The above data set was 2^9 in length. Therefore there are
# coefficients at resolution levels 0, 1, 2, ..., and 8.
#
# The high resolution coefficients are at level 8.
# There should be 256 coefficients at level 8 in index location 0 and 1.
#
length(getpacket(MyWST, level=8, index=0))
# [1] 256
length(getpacket(MyWST, level=8, index=1))
# [1] 256
#
# There should be 4 coefficients at resolution level 2
#
getpacket(MyWST, level=2, index=0)
# [1] -0.92103095  0.70125471  0.07361174 -0.43467375
#
# O.k. Let's insert the packet containing the numbers 19,42,21,32
#
NewMyWST <- putpacket(MyWST, level=2, index=0, packet=c(19,42,31,32))
#
# Let's check that it put the numbers in correctly by reaccessing that
# packet...
#
getpacket(NewMyWST, level=2, index=0)
# [1] 19 42 31 32
#
# Yep. It inserted the packet correctly.
}
\keyword{manip}
\author{G P Nason}
