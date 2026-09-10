\name{accessC.wp}
\alias{accessC.wp}
\title{Warning function when trying to access smooths from wavelet packet object (wp).}
\description{
   There are no real smooths to access in a \code{\link{wp}} wavelet packet object. This function returns an error message. To obtain coefficients from a wavelet packet object you should use the \code{\link{getpacket}} collection of functions.
}
\usage{
\method{accessC}{wp}(wp, \dots)
}
\arguments{
 \item{wp}{Wavelet packet object.}
\item{\dots}{any other arguments}
}
\value{An error message!}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994}
\seealso{ \code{\link{getpacket}}
}
\keyword{manip}
\author{G P Nason}
