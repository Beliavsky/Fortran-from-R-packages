\name{accessD.wpst}
\alias{accessD.wpst}
\title{Get coefficients from a non-decimated wavelet packet object (wpst) in time order.}
\description{
The coefficients from a non-decimated wavelet packet object, \code{\link{wpst}}, are stored in a particular order in the wpst component of the wpstobj object. This function extracts all the coefficients corresponding to a particular wavelet packet in time order.
}
\usage{
\method{accessD}{wpst}(wpst, level, index, \dots)
}
\arguments{
\item{wpst}{Non-decimated wavelet packet object from which you wish to extract time-ordered coefficients.}
\item{level}{The resolution level that you wish to extract. This can range from zero (the coarsest coefficients) to nlevelsWT-1(wstobj) which are the finest scale coefficients.}
\item{index}{The wavelet packet index that you require (sequency ordering). This can range from 0 (father wavelet coeffcients) to \code{2^(nlevelsWT - level) - 1}, i.e. the maximum is dependent on the resolution level.}
\item{\dots}{any other arguments}
}
\details{The \code{\link{wpst}} function performs a non-decimated wavelet packet transform. This function extracts the coefficients at a particular resolution level specified by level in time order. 

It is possible to extract the individual packets (before interweaving, i.e. the direct result of multiple applications of the packet operators) by using the \code{\link{getpacket.wpst}} function.}

\references{
Nason, G.P., Sapatinas, T. and Sawczenko, A. Statistical modelling using undecimated wavelet transforms. 
}
\seealso{
\code{\link{wpst}}, \code{wpst.object}, \code{\link{accessD}}, \code{\link{getpacket.wpst}}
}
\examples{
#
# Get the 4th level of coefficients from a decomposition
#
dat <- rnorm(128)
accessD(wpst(dat), level=4, index=3)
}
\keyword{manip}
\author{G P Nason}
