\name{getpacket.wp}
\alias{getpacket.wp}
\title{Get packet of coefficients from a wavelet packet object (wp).}
\description{
This function extracts and returns a packet of coefficients from a wavelet packet (\code{\link{wp}}) object.
}
\usage{
\method{getpacket}{wp}(wp, level, index, \dots )
}
\arguments{
\item{wp}{Wavelet packet object from which you wish to extract the packet from.} 
\item{level}{The resolution level of the coefficients that you wish to extract.} 
\item{index}{The index number within the resolution level of the packet of coefficients that you wish to extract.} 
\item{\dots}{any other arguments}
}
\details{
The \code{\link{wp}} produces a wavelet packet object. The coefficients in this structure can be organised into a binary tree with each node in the tree containing a packet of coefficients. 

Each packet of coefficients is obtained by chaining together the effect of the \emph{two packet operators} DG and DH: these are the high and low pass quadrature mirror filters of the Mallat pyramid algorithm scheme followed by decimation (see Mallat~(1989b)). 

Starting with data \eqn{c^J} at resolution level J containing
\eqn{2^J} data points the wavelet packet algorithm operates as follows.
First DG and DH are applied to \eqn{c^J} producing
\eqn{d^{J-1}} and \eqn{c^{J-1}} respectively.
Each of these sets of coefficients is of length one half of the original data:
i.e. \eqn{2^{J-1}}.
Each of these sets of coefficients is a set of
\emph{wavelet packet coefficients}.
The algorithm then applies both DG and DH to both
\eqn{d^{J-1}} and \eqn{c^{J-1}} to form a four sets of coefficients at level
J-2. Both operators are used again on the four sets to produce 8 sets, then again on the 8 sets to form 16 sets and so on. At level j=J,...,0 there are
\eqn{2^{J-j}} packets of coefficients each containing \eqn{2^j} coefficients. 

This function enables whole packets of coefficients to be extracted at any resolution level. The index argument chooses a particular packet within each level and thus ranges from 0 (which always refer to the father wavelet coefficients), 1 (which always refer to the mother wavelet coefficients) up to \eqn{2^{J-j}}. 

}
\value{
A vector containing the packet of wavelet packet coefficients that you wished to extract. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{wp}}, \code{\link{putpacket.wp}}, \code{\link{basisplot.wp}}, \code{\link{draw.wp}}, \code{\link{InvBasis.wp}}, \code{\link{MaNoVe.wp}}, \code{nlevelsWT.wp}, \code{\link{plot.wp}}. \code{\link{threshold.wp}}. 
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
#[1] 256
length(getpacket(MyWP, level=8, index=1))
#[1] 256
#
# The next command shows that there are only two packets at level 8
#
\dontrun{getpacket(MyWP, level=8, index=2)}
#Index was too high, maximum for this level is  1 
#Error in getpacket.wp(MyWP, level = 8, index = 2): Error occured
#Dumped
#
# There should be 4 coefficients at resolution level 2
#
# The father wavelet coefficients are (index=0)
getpacket(MyWP, level=2, index=0)
#[1] -0.9736576  0.5579501  0.3100629 -0.3834068
#
# The mother wavelet coefficients are (index=1)
#
#[1]  0.72871405  0.04356728 -0.43175307  1.77291483
#
# There will be 127 packets at this level.
#
}
\keyword{manip}
\author{G P Nason}
