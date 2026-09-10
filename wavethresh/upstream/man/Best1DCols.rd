\name{Best1DCols}
\alias{Best1DCols}
\title{Extract the best (one-dimensional) nondecimated WP packets}
\description{
This function takes the whole set of nondecimated wavelet packets and
selects those packets that correlate best with the "response" groups.
The idea is to reduce the large dimensionality (number of packets) into
something more manageable which can then be fed into a proper discriminator.
}
\usage{
Best1DCols(w2d, mincor=    0.69999999999999996)
}
\arguments{
\item{w2d}{An object that gets returned from a call to the
	\code{\link{wpst2discr}} function which turns a wpst class
	object into a regular multivariate matrix}
\item{mincor}{The threshold above which variables (packets) 
	get included into the
	final mix if their correlation with the groups variable
	 is higher than this value.}
}
\details{This function is not intended for direct user use.
	In this function, the w2d object contains a matrix where each
	column contains the coefficients of a single packet from a
	non-decimated wavelet packet transform. The number of rows of the
	matrix is the same as the original time series and hence each
	column can be correlated with a separate groups variable that
	contains the group membership of a separate variable which changes
	over time. Those packet columns that have correlation greater
	than the \code{mincor} value are extracted and returned
	in the \code{BasisMatrix} item of the returned list. 
}
\value{A list with the following components:
\item{nlevelsWT}{The number of levels of the nondecimated wavelet packet
	encapsulator, w2d}
\item{BasisMatrix}{The highest correlating packets, sorted according to
	decreasing correlation}
\item{level}{The levels corresponding to the selected packets}
\item{pkt}{The packet indices corresponding to the selected packets}
\item{basiscoef}{The sorted correlations}
\item{groups}{The groups time series}
}
\seealso{\code{\link{makewpstDO}},\code{\link{wpst2discr}}}
\author{G P Nason}
\keyword{ts}
\keyword{multivariate}
