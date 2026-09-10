\name{wpst}
\alias{wpst}
\title{Non-decimated wavelet packet transform. }
\description{
This function computes the non-decimated wavelet packet transform as described by Nason, Sapatinas and Sawczenko, 1998. The non-decimated wavelet packet transform (NWPT) contains all possible shifted versions of the wavelet packet transform.
}
\usage{
wpst(data, filter.number=10, family="DaubLeAsymm", FinishLevel)
}
\arguments{
\item{data}{A vector containing the data you wish to decompose. The length of this vector must be a power of 2.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{FinishLevel}{At which level to stop decomposing. The full decomposition
	decomposes to level 0, but you could stop earlier.} 
}
\details{
This function computes the packet-ordered non-decimated wavelet packet transform of data as described by Nason, Sapatinas and Sawczenko, 1998. It assumes periodic boundary conditions. The order of computation of the NWPT is
\eqn{O(n^2)} if n is the number of input data points. 

Packets can be extracted from the \code{wpst.object} produced by this function using the \code{\link{getpacket.wpst}} function. Whole resolution levels of non-decimated wavelet packet coefficients in time order can be obtained by using the \code{\link{accessD.wpst}} function. 
}
\value{
An object of class \code{\link{wpst}} containing the discrete packet-ordered non-decimated wavelet packet coefficients. 
}
\section{RELEASE}{Version 3.8.8 Copyright Guy Nason 1997 }
\seealso{
\code{\link{accessD}}, \code{\link{accessD.wpst}}, \code{\link{filter.select}}, \code{\link{getpacket}}, \code{\link{getpacket.wpst}},
\code{\link{makewpstDO}}
}
\examples{
v <- rnorm(128)
vwpst <- wpst(v)
}
\keyword{math}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
