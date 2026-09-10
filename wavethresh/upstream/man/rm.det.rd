\name{rm.det}
\alias{rm.det}
\title{Set coarse levels of a wavelets on the interval transform object to zero}
\usage{
rm.det(wd.int.obj)
}
\arguments{
\item{wd.int.obj}{the object whose coarse levels you wish to set to zero}
}
\description{
Set the wavelet coefficients of certain coarse levels for a "wavelets on
the interval" object equal to zero. The operation of this function
is somewhat similar to the \code{\link{nullevels}} function, but for
objects associated with the "wavelets on the interval code".
}
\details{
The "wavelets on the interval" code is contained within the \code{\link{wd}}
function. All levels coarser than (but not including) the
\code{wd.int.obj$current.scale} are set to zero.
}
\value{
A \code{\link{wd.object}} of \code{type="interval"} containing the modified
input object with certain coarse levels set to zero.
}
\seealso{\code{\link{nullevels}}, \code{\link{wd}}}
\author{Piotr Fryzlewicz}
\keyword{manip}
