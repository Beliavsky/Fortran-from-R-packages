\name{drawwp.default}
\alias{drawwp.default}
\title{Subsidiary routine that actually computes wavelet packet values}
\usage{
drawwp.default(level, index, filter.number = 10, family = "DaubLeAsymm",
    resolution = 64 * 2^level)
}
\arguments{
\item{level}{The resolution level of the packet you want}
\item{index}{The packet index of the packet you want}
\item{filter.number}{The type of wavelet you want, see \code{\link{filter.select}}}
\item{family}{The family of wavelet you want, see \code{\link{filter.select}}}
\item{resolution}{The number of ordinates at which you want the wavelet packet}
}
\description{
Function computes the values of a given wavelet packet on a discrete
grid.
}
\details{
Function works by computing a wavelet packet transform of a zero vector.
Then inserting a single one somewhere in the desired packet, and then
inverts the transform.
}
\value{
A vector containing the "y" values of the required wavelet packet. 
}
\seealso{\code{\link{draw.wp}},\code{\link{InvBasis}},
\code{\link{nlevelsWT}},
\code{\link{putpacket}},
\code{\link{wp}}}
\author{G P Nason}
\keyword{dplot}
