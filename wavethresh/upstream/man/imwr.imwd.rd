\name{imwr.imwd}
\alias{imwr.imwd}
\title{Inverse two-dimensional discrete wavelet transform.}
\description{
This functions performs the reconstruction stage of Mallat's pyramid algorithm (i.e. the inverse discrete wavelet transform) for images. 
}
\usage{
\method{imwr}{imwd}(imwd, bc=imwd$bc, verbose=FALSE, \dots) 
}
\arguments{
\item{imwd}{An object of class `\code{\link{imwd}}'. This type of object is returned by `\code{\link{imwd}}'.}
\item{bc}{This argument specifies the boundary handling, it is best left to be the boundary handling specified by that in the supplied imwd (as is the default).}
\item{verbose}{If this argument is true then informative messages are printed detailing the computations to be performed}
\item{\dots}{any other arguments}
}
\details{
Details of the algorithm are to be found in Mallat (1989). Similarly to the decomposition function, \code{\link{imwd}} the inverse algorithm works by applying many 1D reconstruction algorithms to the coefficients. The filters in these 1D reconstructions are incorporated in the supplied \code{\link{imwd.object}} and originally created by the \code{\link{filter.select}} function in WaveThresh3.
 
This function is a method for the generic function \code{\link{imwr}} for class \code{\link{imwd.object}}. It can be invoked by calling \code{\link{imwr}} for an object of the appropriate class, or directly by calling imwr.imwd regardless of the class of the object. 
}
\value{
A matrix, of dimension determined by the original data set supplied to the initial decomposition (more precisely, determined by the \code{\link{nlevelsWT}} component of the \code{\link{imwd.object}}). This matrix is the highest resolution level of the reconstruction. If a \code{\link{imwd}} two-dimensional wavelet transform is followed immediately by a \code{\link{imwr}} inverse two-dimensional wavelet transform then the returned matrix will be exactly the same as the original image. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{imwr}}. 
}
\examples{
#
# Do a decomposition, then exact reconstruction
# Look at the error
#
test.image <- matrix(rnorm(32*32), nrow=32)
#
# Test image is just some sort of  square matrix whose side length
# is a power of two.
#
max( abs(imwr(imwd(test.image)) - test.image))
# [1] 1.014611e-11
}
\keyword{nonlinear}
\keyword{smooth}
\author{G P Nason}
