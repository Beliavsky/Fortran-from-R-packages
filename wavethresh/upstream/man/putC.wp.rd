\name{putC.wp}
\alias{putC.wp}
\title{Warning function when trying to insert father wavelet coefficients into wavelet packet object (wp).}
\description{
There are no real smooths to insert in a \code{\link{wp}} wavelet packet object. This function returns an error message. To insert coefficients into a wavelet packet object you should use the \code{\link{putpacket}} collection of functions. 
}
\usage{
\method{putC}{wp}(wp, \dots)
}
\arguments{
\item{wp}{Wavelet packet object.} 
\item{\dots}{any other arguments}
}
\value{
An error message! 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{putpacket}}, \code{\link{putpacket.wp}}. 
}
\keyword{error}
\author{G P Nason}
