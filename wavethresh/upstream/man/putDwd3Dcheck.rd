\name{putDwd3Dcheck}
\alias{putDwd3Dcheck}
\title{Check argument list for putD.wd3D}
\description{
This function checks the argument list for \code{\link{putD.wd3D}} and is not meant to be directly called by any user.
}
\usage{
putDwd3Dcheck(lti, dima, block, nlx)
}
\arguments{
\item{lti}{At which level of the \code{\link{wd3D.object}} you wish to insert a block of coefficients.}
\item{dima}{A vector, of length 3, which specifies the dimension of the block to insert.} 
\item{block}{A character string which specifies which block is being inserted (one of GGG, GGH, GHG, GHH, HGG, HGH, HHG, or HHH).}
\item{nlx}{The number of levels in the \code{\link{wd3D.object}} that you wish to insert the coefficients into (can be obtained using the \code{\link{nlevelsWT}} function). }
}
\details{
This function merely checks that the dimensions and sizes of the array to be inserted into a \code{\link{wd3D.object}} using the \code{\link{putD.wd3D}} function are correct. 
}
\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{accessD}}, \code{\link{putD}}, \code{\link{accessD.wd3D}}, \code{\link{print.wd3D}}, \code{\link{putD}}, \code{\link{summary.wd3D}}, \code{\link{threshold.wd3D}}, \code{\link{wd3D}}, \code{\link{wd3D.object}}, \code{\link{wr3D}}.
 }
\examples{
#
# Not intended to be used by the user!
#
}
\keyword{manip}
\author{G P Nason}
