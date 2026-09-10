\name{accessD.wp}
\alias{accessD.wp}
\title{Obtain whole resolution level of wavelet packet coefficients from a wavelet packet object (wp).}
\description{
Get a whole resolution level's worth of coefficients from a \code{\link{wp}} wavelet packet object. To obtain packets of coefficients from a wavelet packet object you should use the \code{\link{getpacket}} collection of functions.
}
\usage{
\method{accessD}{wp}(wp, level, \dots)
}
\arguments{
\item{wp}{Wavelet packet object}.
\item{level}{the resolution level that you wish to extract.}
\item{\dots}{any other arguments}
}
\details{
The wavelet packet coefficients are actually stored in a straightforward manner in a matrix component of a \code{\link{wp}} object so it would not be too difficult to extract whole resolution levels yourself. However, this routine makes it easier to do.
}
\value{
A vector containing the coefficients that you wanted to extract. 
}
\section{RELEASE}{
Version 3.5.3 Copyright Guy Nason 1994 
}
\seealso{
\code{\link{accessD}}, \code{\link{getpacket}}
}
\keyword{manip}
\author{G P Nason}
