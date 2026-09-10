\name{guyrot}
\alias{guyrot}
\title{Cyclically rotate elements of a vector}
\usage{
guyrot(v, n)
}
\arguments{
\item{v}{Vector whose elements you wish to rotate}
\item{n}{Integer determining the amount to rotate, can be negative}
}
\description{
This function shifts (or rotates) the elements of the input vector in a cyclic fashion (end periodicity is used).
}
\details{
A very simple function which cyclically shifts the elements of a vector. Not necessarily intended as a top level user function but it is a useful little function.
}
\value{
A vector containing the shifted or rotated coefficients.
}
\seealso{\code{\link{wpst2discr}},
\code{\link{wpstCLASS}}}
\examples{
#
# Start off with an example vector
#
v <- c(1,2,3,4,5,6)
#
# Rotate it one element to the right, rightmost element gets rotated round
# to be first element.
#
guyrot(v,1)
# [1] 6 1 2 3 4 5
#
# Rotate v two spaces to the left, leftmost two elements get rotated around
# to be new last elements
guyrot(v, -2)
#
# [1] 3 4 5 6 1 2
#
#
# Now issue a larger rotation, e.g. 19!
#
guyrot(v,19)
# [1] 6 1 2 3 4 5
#
# Its just the same as rotating by 1 since the input vector is of length 6
# and so rotating by 19 is the same as rotating by 6,6,6, and then 1!
#
}
\author{G P Nason}
\keyword{math}
