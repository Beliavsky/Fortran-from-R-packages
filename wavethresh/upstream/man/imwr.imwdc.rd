\name{imwr.imwdc}
\alias{imwr.imwdc}
\title{Inverse two-dimensional discrete wavelet transform. }
\description{
Inverse two-dimensional discrete wavelet transform. 
}
\usage{
\method{imwr}{imwdc}(imwd, verbose=FALSE, \dots)
}
\arguments{
\item{imwd}{An object of class \code{imwdc}. This type of object is returned by \code{\link{threshold.imwd}} and is a \code{\link{compress.imwd}} compressed version of an \code{\link{imwd}} object.}
\item{verbose}{If this argument is true then informative messages are printed detailing the computations to be performed}
\item{\dots}{other arguments to supply to the \code{\link{imwr}} function which is called after uncompressing the imwdc object.}} 

\details{
This function merely uncompresses the supplied \code{\link{imwdc.object}} and passes the resultant \code{\link{imwd}} object to the \code{\link{imwr.imwd}} function.
 
This function is a method for the generic function \code{\link{imwr}} for class \code{\link{imwdc.object}}. It can be invoked by calling \code{\link{imwr}} for an object of the appropriate class, or directly by calling imwr.imwdc regardless of the class of the object. 
}
\value{
A matrix, of dimension determined by the original data set supplied to the initial decomposition (more precisely, determined by the \code{\link{nlevelsWT}} component of the \code{\link{imwdc.object}}). This matrix is the highest resolution level of the reconstruction. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{compress.imwd}}, \code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{imwr}}. 
}
\examples{
#
# Do a decomposition, thresholding, then exact reconstruction
# Look at the error
#

test.image <- matrix(rnorm(32*32), nrow=32)

# Test image is just some sort of  square matrix whose side length
# is a power of two.
#
max( abs(imwr(threshold(imwd(test.image))) - test.image))
# [1] 62.34
#
# The answer is not zero (see contrasting examples in the help page for
# imwr.imwd because we have thresholded the
# 2D wavelet transform here).
}
\keyword{nonlinear}
\keyword{smooth}
\author{G P Nason}
