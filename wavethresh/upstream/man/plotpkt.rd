\name{plotpkt}
\alias{plotpkt}
\title{Sets up a high level plot ready to show the time-frequency plane
	and wavelet packet basis slots}
\usage{
plotpkt(J)
}
\arguments{
\item{J}{The number of resolution levels associated with the wavelet
	packet object you want to depict}
}
\description{
Sets up a high level plot ready to add wavelet packet slots using, 
e.g. \code{\link{addpkt}}. This function is used by several routines
to begin plotting graphical representations of the time-frequency plane
and spaces for packets.
}
\details{
Description says all
}
\value{
Nothing of interest
}
\seealso{\code{\link{addpkt}}, \code{\link{basisplot}}, \code{\link{basisplot.BP}}, \code{\link{basisplot.wp}}, \code{\link{plot.nvwp}}}
\author{G P Nason}
\keyword{dplot}
