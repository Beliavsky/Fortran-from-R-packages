\name{getpacket.wst}
\alias{getpacket.wst}
\title{Get packet of coefficients from a packet ordered non-decimated wavelet object (wst).}
\description{
This function extracts and returns a packet of coefficients from a packet-ordered non-decimated wavelet object (\code{\link{wst}}) object. The \code{\link{wst}} objects are computed by the \code{\link{wst}} function amongst others. 
}
\usage{
\method{getpacket}{wst}(wst, level, index, type="D", aspect, \dots)
}
\arguments{
\item{wst}{Packet-ordered non-decimated wavelet object from which you wish to extract the packet from.} 
\item{level}{The resolution level of the coefficients that you wish to extract.} 
\item{index}{The index number within the resolution level of the packet of coefficients that you wish to extract.} 
\item{type}{This argument must be either "\code{C}" or "\code{D}". If the argument is "\code{C}" then non-decimated father wavelet coefficients corresponding to the packet that you want are returned. If the argument is "\code{D}" then non-decimated mother wavelet coefficients are returned. }
\item{aspect}{Function applied to the coefficients before return.
	This is suppled as a character string which gets converted to
	a function to apply. For example, "Mod" for complex-valued
	coefficients returns the absolute values.}
\item{\dots}{Other arguments}
}
\details{
The \code{\link{wst}} function produces a packet-ordered non-decimated wavelet object: \code{\link{wst}}. The coefficients in this structure can be organised into a binary tree with each node in the tree containing a packet of coefficients. 

Each packet is obtained by repeated application of the usual DG quadrature mirror filter with both even and odd dyadic decimation. See the detailed description given in Nason and Silverman, 1995. 

This function enables whole packets of coefficients to be extracted at any resolution level. The index argument chooses a particular packet within each level and thus ranges from 0 to \eqn{2^{J-j}} for j=0,..., J-1. Each packet corresponds to the wavelet coefficients with respect to different origins. 

Note that both mother and father wavelet coefficient at different shifts are available by using the type argument. 
}
\value{
A vector containing the packet of packet-ordered non-decimated wavelet coefficients that you wished to extract. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, 
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
#[1] 256
length(getpacket(MyWST, level=8, index=1))
#[1] 256
#
# There are also 256 FATHER wavelet coefficients at each of these two indices
# (origins)
#
length(getpacket(MyWST, level=8, index=0, type="C"))
#[1] 256
length(getpacket(MyWST, level=8, index=1, type="C"))
#[1] 256
#
# There should be 4 coefficients at resolution level 2
#
getpacket(MyWST, level=2, index=0)
#[1] -0.92103095  0.70125471  0.07361174 -0.43467375
#
# Here are the equivalent father wavelet coefficients
#
getpacket(MyWST, level=2, index=0, type="C")
#[1] -1.8233506 -0.2550734  1.9613138  1.2391913
}
\keyword{manip}
\author{G P Nason}
