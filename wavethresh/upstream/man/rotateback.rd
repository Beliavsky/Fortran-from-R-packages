\name{rotateback}
\alias{rotateback}
\title{Cyclically shift a vector one place to the right}
\usage{
rotateback(v)
}
\arguments{
\item{v}{The vector to shift}
}
\description{
Cyclically shifts the elements of a vector one place to the right. The right-most element becomes the first element.
}
\details{
Subsidiary function used by the \code{\link{av.basis}}
function which is the R function component of the
\code{\link{AvBasis.wst}} function.
}
\value{
The rotated vector
}
\examples{
#
# Here is a test vector
#
v <- 1:10
#
# Apply this function
#
rotateback(v)
#[1] 10  1  2  3  4  5  6  7  8  9
#
# A silly little function really!
}
\author{G P Nason}
\keyword{math}
