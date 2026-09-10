\name{av.basis}
\alias{av.basis}
\title{Perform basis averaging for wst class object}
\usage{
av.basis(wst, level, ix1, ix2, filter)
}
\arguments{
\item{wst}{The \code{\link{wst.object}} that you wish to basis average}
\item{level}{The resolution level the function is currently operating at}
\item{ix1}{Which "left" packet in the level you are accessing}
\item{ix2}{Which "right" packet}
\item{filter}{The wavelet filter details, see \code{\link{filter.select}}}
}
\description{
\bold{Note:} that this function is not for direct user use. This function
is a helper routine for the \code{\link{AvBasis.wst}} function which is
the one that should be used by users.

This function works by recursion, essentially it merges the current
levels C coefficients from one packet shift with its associated
D coefficients, does the same for the other packet shift and then
averages the two reconstructions to provide the C coefficients for
the next level up.
}
\details{
Description says all, see help page for \code{\link{AvBasis.wst}}.
}
\value{
Returns the average basis reconstruction of a \code{\link{wst.object}}.
}
\seealso{\code{\link{AvBasis}}, \code{\link{AvBasis.wst}},
	\code{\link{conbar}}, \code{\link{rotateback}},
	\code{\link{getpacket}}}
\author{G P Nason}
\keyword{manip}
