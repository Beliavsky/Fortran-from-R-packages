\name{compress.default}
\alias{compress.default}
\title{Do "zero" run-length encoding compression of a vector of numbers. }
\description{
Efficiently compress a vector containing many zeroes. 
}
\usage{
\method{compress}{default}(v, verbose=FALSE,\dots)

}
\arguments{
\item{v}{The vector that you wish to compress. This compression function is efficient at compressing vectors with many zeroes, but is not a \emph{general} compression routine.}
\item{verbose}{If\code{TRUE} then this routine prints out the degree of compression achieved. }
\item{\dots}{any other arguments}
}
\details{
Images are large objects. Thresholded 2d wavelet objects (\code{\link{imwd}}) are also large, but many of their elements are zero. compress.default takes a vector, decides whether compression is necessary and if it is makes an object of class \code{compressed} containing the nonzero elements and their position in the original vector. 

The decision whether to compress the vector or not depends on two things, first the number of non-zero elements in the vector (r, say), and second the length of the vector (n, say). Since the position and value of the non-zero elements is stored we will need to store 2r values for the non-zero elements. So compression takes place if \code{2r < n}. 

This function is the default method for the generic function \code{\link{compress}}. It can be invoked by calling compress for an object of the appropriate class, or directly by calling compress.default regardless of the class of the object. 

}
\value{
An object of class compressed if \code{compression} took place, otherwise a an object of class \code{uncompressed}. 
}
\note{
Sometimes the compressed object can be larger than the original. This usually only happens for small objects, so it doesn't really matter. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{compress}},  \code{\link{imwd}}, \code{\link{threshold.imwd}},
\code{\link{uncompress}}
}
\examples{
#
# Compress a vector with lots of zeroes
#
compress(c(rep(0,100),99))
#$position:
#[1] 101
#
#$values:
#[1] 99
#
#$original.length:
#[1] 101
#
#attr(, "class"):
#[1] "compressed"
#
# Try to compress a vector with not many zeroes
#
compress(1:10)
#$vector:
#[1]  1 2  3  4  5  6  7  8  9 10
#
#attr(, "class"):
#[1] "uncompressed"
#
#
}
\keyword{manip}
\author{G P Nason}
