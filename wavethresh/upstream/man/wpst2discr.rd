\name{wpst2discr}
\alias{wpst2discr}
\title{Reshape/reformat packet coefficients into a multivariate data set}
\description{
The packet coefficients of a nondecimated wavelet packet object are stored
internally in an efficient form. This function takes the nondecimated
wavelet packets and stores them as a matrix (multivariate data set).
Each column in the returned matrix corresponds to an individual packet,
each row corresponds to a time index in the original packet or time series.
}
\usage{
wpst2discr(wpstobj, groups)
}
\arguments{
\item{wpstobj}{A wpst class object, output from \code{\link{wpst}} say}
\item{groups}{A time series containing the group membership at each time point}
}
\details{
Description says it all
}
\value{An object of class w2d which is a list containing the following
items:
\item{m}{The matrix containing columns of packet information.}
\item{groups}{Passes through the \code{group} argument from input.}
\item{level}{Each column corresponds to a packet, this vector contains the
	information on which resolution level each packet comes from}
\item{pktix}{Like for \code{level} but for packet indices}
\item{nlevelsWT}{The number of resolution levels in total, from the wpst object}
}
\seealso{\code{\link{makewpstDO}}, \code{\link{wpst}}}
\author{G P Nason}
\keyword{multivariate}
\keyword{ts}
