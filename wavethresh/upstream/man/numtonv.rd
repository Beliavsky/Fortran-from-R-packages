\name{numtonv}
\alias{numtonv}
\title{Convert an index number into a node vector object.}
\description{
Convert an index number into a \code{node vector} object. 
}
\usage{
numtonv(number, nlevels)
}
\arguments{
\item{number}{The index number of a particular basis within a wavelet object.} 
\item{nlevels}{The number of levels that the wavelet object has (can often be discovered using the \code{\link{nlevels}} function). }
}
\details{
A basis within a (e.g. non-decimated) wavelet object (such as a \code{\link{wst.object}}) is represented in WaveThresh by a \code{nv} or node vector. 

A packet-ordered non-decimated wavelet transform object \code{\link{wst}} for short) which is the transform of a vector of length \code{n} contains \code{n} bases. Each basis can be indexed from 0 to \code{(n-1)} . 

A \code{\link{wst.object}} is simply a fully populated binary tree. There are nlevels levels in the tree with a split at each level. The root of the tree is at level 0, there are two branches at level 1, four at level 2, eight at level 3 and so on. A path through the tree can be constructed by starting at the root and choosing "left" or "right" at each possible branch. For certain data situations this path is constructed using minimum entropy algorithms (for examples \code{\link{MaNoVe}}). This function (numtonv takes the numerical representation of a path and converts it into a \code{node.vector} form suitable for passing to \code{\link{InvBasis}} to invert the representation according to a basis specicified by number.
 
The least significant digit in number corresponds to deciding on the left/right decision at the fine leaves of the tree (high-frequency structure) and the
most significant digit in number corresponds to deciding on the left/right
decision at the root.
Therefore gradually incrementing number from 0 to
\code{2^{nlevels}-1} steps through all possible bases in the
\code{\link{wst}} object ranging from all decisions being made "left" to all decisions being made "right".
The "number" dividied by \code{2^{nlevels}} corresponds exactly to the binary number epsilon in Nason and Silverman (1995). 
}
\value{
An object of class \code{nv} (node vector). This contains information about a path through a wavelet object (a basis in a wavelet object). 
}
\section{RELEASE}{Version 3.6.0 Copyright Guy Nason 1995 }
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{MaNoVe}}, \code{\link{nv.object}}, \code{\link{InvBasis}}, \code{\link{nlevels}}.
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
#
# Make it noisy
#
ynoise <- test.data + rnorm(512, sd=0.1)
#
# Do packet ordered non-decimated wavelet transform
#
ynwst <- wst(ynoise)
#
# Now threshold the coefficients
#
ynwstT <- threshold(ynwst)
#
# Select basis number 9 (why not?)
#
NodeVector9 <- numtonv(9, nlevelsWT(ynwstT))
#
# Let's print it out to see what it looks like
# (nb, if you're repeating this examples, the basis might be different
# as you may have generated different pseudo random noise to me)
#
NodeVector9
# Level :  8  Action is  R (getpacket Index:  1 )
# Level :  7  Action is  L (getpacket Index:  2 )
# Level :  6  Action is  L (getpacket Index:  4 )
# Level :  5  Action is  R (getpacket Index:  9 )
# Level :  4  Action is  L (getpacket Index:  18 )
# Level :  3  Action is  L (getpacket Index:  36 )
# Level :  2  Action is  L (getpacket Index:  72 )
# Level :  1  Action is  L (getpacket Index:  144 )
# Level :  0  Action is  L (getpacket Index:  288 )
# There are  9  reconstruction steps
#
# The print-out describes the tree through ynwstT that corresponds to
# basis 9.
#
# The NodeVector9 and ynwstT objects could now be supplied to
# InvBasis.wst for inverting ynwstT according
# to the NodeVector9 or basis number 9.
}
\keyword{algebra}
\author{G P Nason}
