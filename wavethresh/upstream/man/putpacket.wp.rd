\name{putpacket.wp}
\alias{putpacket.wp}
\title{Inserts a packet of coefficients into a wavelet packet object (wp).}
\description{
This function inserts a packet of coefficients into a wavelet packet (\code{\link{wp}}) object.
}
\usage{
\method{putpacket}{wp}(wp, level, index, packet , \dots)
}
\arguments{
\item{wp}{Wavelet packet object into which you wish to put the packet.} 
\item{level}{The resolution level of the coefficients that you wish to insert.}
\item{index}{The index number within the resolution level of the packet of coefficients that you wish to insert.}
\item{packet}{a vector of coefficients which is the packet you wish to insert.} \item{\dots}{any other arguments}
}
\details{
The coefficients in this structure can be organised into a binary tree with each node in the tree containing a packet of coefficients. 

Each packet of coefficients is obtained by chaining together the effect of the two \emph{packet operators} DG and DH: these are the high and low pass quadrature mirror filters of the Mallat pyramid algorithm scheme followed by decimation (see Mallat (1989b)). 

Starting with data \eqn{c^J} at resolution level J containing
\eqn{2^J} data points the wavelet packet algorithm operates as follows.
First DG and DH are applied to
\eqn{c^J} producing \eqn{d^{J-1}} and \eqn{c^{J-1}} respectively.
Each of these sets of coefficients is of length one half of the original data: i.e. \eqn{2^{J-1}}. Each of these sets of coefficients is a set of \emph{wavelet packet coefficients}. The algorithm then applies both DG and DH to both
\eqn{d^{J-1}} and \eqn{c^{J-1}} to form a four sets of coefficients at
level J-2. Both operators are used again on the four sets to produce 8 sets, then again on the 8 sets to form 16 sets and so on.
At level j=J,...,0 there are \eqn{2^{J-j}} packets of coefficients each
containing \eqn{2^j} coefficients. 

This function enables whole packets of coefficients to be inserted at any resolution level. The \code{index} argument chooses a particular packet within each level and thus ranges from 0 (which always refer to the father wavelet coefficients), 1 (which always refer to the mother wavelet coefficients) up to \eqn{2^{J-j}}.

}
\value{
An object of class \code{\link{wp.object}} which is the same as the input \code{\link{wp.object}} except it now has a modified packet of coefficients. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{wp}}, \code{\link{getpacket.wp}}. \code{\link{putpacket}}. 
}
\examples{
#
# Take the wavelet packet transform of some random data
#
MyWP <- wp(rnorm(1:512))
#
# The above data set was 2^9 in length. Therefore there are
# coefficients at resolution levels 0, 1, 2, ..., and 8.
#
# The high resolution coefficients are at level 8.
# There should be 256 DG coefficients and 256 DH coefficients
#
length(getpacket(MyWP, level=8, index=0))
# [1] 256
length(getpacket(MyWP, level=8, index=1))
# [1] 256
#
# The next command shows that there are only two packets at level 8
#
#getpacket(MyWP, level=8, index=2)
# Index was too high, maximum for this level is  1 
# Error in getpacket.wp(MyWP, level = 8, index = 2): Error occured
# Dumped
#
# There should be 4 coefficients at resolution level 2
#
# The father wavelet coefficients are (index=0)
getpacket(MyWP, level=2, index=0)
# [1] -0.9736576  0.5579501  0.3100629 -0.3834068
#
# The mother wavelet coefficients are (index=1)
#
getpacket(MyWP, level=2, index=1)
# [1]  0.72871405  0.04356728 -0.43175307  1.77291483
#
# Well, that exercised the getpacket.wp
# function. Now that we know that level 2 coefficients have 4 coefficients
# let's insert some into the MyWP object.
#
MyWP <- putpacket(MyWP, level=2, index=0, packet=c(21,32,67,89))
#
# O.k. that was painless. Now let's check that the correct coefficients
# were inserted.
#
getpacket(MyWP, level=2, index=0)
#[1] 21 32 67 89
#
# Yep. The correct coefficients were inserted.
}
\keyword{manip}
\author{G P Nason}
