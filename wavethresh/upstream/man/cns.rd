\name{cns}
\alias{cns}
\title{Create new zeroed spectrum.}
\description{
Part of a two-stage function suite designed to simulate locally stationary wavelet processes in conjunction with the LSWsim function. 
}
\usage{
cns(n, filter.number=1, family="DaubExPhase")
}
\arguments{
\item{n}{The length of the simulated process that you want to produce. Must be a power of two (for this software).}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments. }
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
}
\details{
This simple routine merely computes the time-ordered non-decimated wavelet transform of a zero vector of the same length as the eventual simulated series that you wish to produce. 

If you look at this routine you will see that it is extremely simple. First, it checks to see whether the n that you supplied is a power of two. If it is then it creates a zero vector of that length. This is then non-decimated wavelet transformed with the appropriate wavelet. 

The output can then be processed and then finally supplied to LSWsim for process simulation.
}
\value{
An object of class: \code{\link{wd}}, and, in fact, of the non-decimated variety. All wavelet coefficients of this are zero.  
}
\seealso{
\code{\link{LSWsim}}, \code{\link{ewspec}} 
}
\keyword{manip}

\author{G P Nason}
