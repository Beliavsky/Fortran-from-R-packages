\name{wp}
\alias{wp}
\title{Wavelet packet transform. }
\description{
This function computes a wavelet packet transform (computed by the complete binary application of the DH and DG packet operators, as opposed to the Mallat discrete wavelet transform which only recurses on the DH operator [low pass]). 
}
\usage{
wp(data, filter.number=10, family="DaubLeAsymm", verbose=FALSE)
}
\arguments{
\item{data}{A vector containing the data you wish to decompose. The length of this vector must be a power of 2.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{verbose}{if \code{TRUE} then (un)helpful messages are printed during the execution.}
}
\details{
The paper by Nason, Sapatinas and Sawczenko, 1998 details this implementation of the wavelet packet transform. A more thorough reference is Wickerhauser, 1994. 
}
\value{
An object of class \code{\link{wp}} which contains the (decimated) wavelet packet coefficients. 
}
\section{RELEASE}{Version 3.0 Copyright Guy Nason 1994 }
\seealso{
\code{\link{accessC.wp}}, \code{\link{accessD.wp}}, \code{\link{basisplot.wp}}, \code{\link{draw.wp}},\code{\link{drawwp.default}}, \code{\link{filter.select}}, \code{\link{getpacket.wp}}, \code{\link{InvBasis.wp}}, \code{\link{MaNoVe.wp}}, \code{\link{plot.wp}}, \code{\link{print.wp}}, \code{\link{putC.wp}}, \code{\link{putD.wp}}, \code{\link{putpacket.wp}}, \code{\link{summary.wp}}, \code{\link{threshold.wp}}, \code{\link{wp.object}}.
 }
\examples{
v <- rnorm(128)
vwp <- wp(v)

}
\keyword{math}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
