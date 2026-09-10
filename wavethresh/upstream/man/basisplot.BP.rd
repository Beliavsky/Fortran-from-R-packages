\name{basisplot.BP}
\alias{basisplot.BP}
\title{Plot time-frequency plane and basis slots associated with basis object}
\usage{
\method{basisplot}{BP}(x, num=min(10, length(BP$level)), ...)
}
\arguments{
\item{x}{The \code{BP} class object, possibly coming from the BP component
	of the object returned by \code{\link{makewpstDO}} that you wish
	to plot}
\item{num}{The number of packets that you wish to add to the plot}
\item{\dots}{Other arguments}
}
\description{
The \code{x} objects store basis information obtained through the
\code{\link{makewpstDO}} object. This function plots where the basis packets
are on the time frequency plane.
}
\details{
Description says all
}
\value{
Nothing of note
}
\seealso{\code{\link{makewpstDO}},\code{\link{Best1DCols}}}
\examples{
#
# See example in help for \code{\link{makewpstDO}}
#
}
\author{G P Nason}
\keyword{hplot}
