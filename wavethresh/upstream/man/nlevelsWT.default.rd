\name{nlevelsWT.default}
\alias{nlevelsWT.default}
\title{Returns number of levels associated with an object}
\description{
This function returns the number of scale levels associated with either a wavelet type object or an atomic object. 
}
\usage{
\method{nlevelsWT}{default}(object, \dots)
}
\arguments{
\item{object}{An object for which you wish to determine how many levels it has or is associated with.} 
\item{\dots}{any other arguments}
}
\details{
This function first checks to see whether the input object has a component called nlevelsWT. If it does then it returns the value of this component. If it does not then it takes the length of the object and then uses the \code{\link{IsPowerOfTwo}} function to return the power of two which equals the length (if any) or NA if the length of the object is not a power of two. 
}
\value{
The number of resolution (scale) levels associated with the object. 
}
\author{Version 3.6.0 Copyright Guy Nason 1995 }
\seealso{
\code{\link{nlevelsWT}}
}
\examples{
#
# Generate some test data 
#
test.data <- example.1()$y
#
# Now, this vector is 512 elements long. What number of levels would any
# wavelet object be that was associated with this vector?
#
nlevelsWT(test.data)
# [1] 9
#
# I.e. 2^9=512. Let's check by taking the wavelet transform of the
# test data and seeing how many levels it actually has
#
nlevelsWT(wd(test.data))
# [1] 9

}
\keyword{arith}
