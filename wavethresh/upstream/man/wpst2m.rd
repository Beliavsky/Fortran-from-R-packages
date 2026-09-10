\name{wpst2m}
\alias{wpst2m}
\title{
Converts a nondecimated wavelet packet object to a (large) matrix with
packets stored as columns
}
\description{
Takes a nondecimated wavelet packet transform, takes the packets one packet
at a time and stores them in a matrix. The packets are rotated on extraction
and storage in the matrix in an attempt to align them, they are also
optionally transformed by \code{trans}. The rotation is performed
by \code{\link{compgrot}}. 

Note that the coefficients are of some series, not the basis functions
themselves. 
}
\usage{
wpst2m(wpstobj, trans = identity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wpstobj}{
	The nondecimated wavelet packet object to store
}
  \item{trans}{
	The optional transform to apply to the coefficients
}
}
\details{
Description says all
}
\value{
A list, of class w2m, with the following components:
\item{m}{The matrix containing the packets}
\item{level}{A vector containing the levels from where the packets in m
	come from}
\item{pktix}{A vector containing the packet indices from where the packets in m
	come from}
\item{nlevelsWT}{The number of resolution levels from the original wpst object}
}
\author{
G P Nason
}

\seealso{\code{\link{makewpstRO}}, \code{\link{print.w2m}}}
\examples{
#
# Not intended to be directly used by users
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
