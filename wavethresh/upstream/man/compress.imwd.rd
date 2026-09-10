\name{compress.imwd}
\alias{compress.imwd}
\title{Compress a (thresholded) imwd class object by removing zeroes. }
\description{
Compress a (thresholded) \code{imwd} class object by removing zeroes. 
}
\usage{
\method{compress}{imwd}(x, verbose=FALSE, \dots) 
}
\arguments{
\item{x}{Object to compress. Compression only does anything on \code{thresholded} \code{\link{imwd.object}}. }
\item{verbose}{If this is true then report on compression activity.}
\item{\dots}{any other arguments}
}
\details{
Thresholded \code{\link{imwd}} objects are usually very large and contain many zero elements. This function compresses these objects into smaller \code{\link{imwd}} objects by using the \code{\link{compress.default}} function which removing the zeroes. This function is a method for the generic function \code{\link{compress}} for class \code{\link{imwd}} objects. It can be invoked by calling \code{\link{compress}} for an object of the appropriate class, or directly by calling \code{\link{compress.imwd}} regardless of the class of the object
}
\value{
An object of type "\code{imwdc}" representing the compressed imwd object. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{compress}}, \code{\link{compress.default}}, \code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{imwdc.object}}, \code{\link{threshold.imwd}}. 
}
\examples{
#
# The user shouldn't need to use this function directly as the
# \code{\link{threshold.imwd}} function calls it
# automatically.
#
}
\keyword{manip}
\author{G P Nason}
