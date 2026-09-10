\name{lt.to.name}
\alias{lt.to.name}
\title{Convert desired level and orientation into code used by imwd}
\usage{
lt.to.name(level, type)
}
\arguments{
\item{level}{Resolution level of coefficients that you want to extract or
	manipulate.}
\item{type}{One of CC, CD, DC or DD indicating smoothed, horizontal,
	vertical or diagonal coefficients}
}
\description{
Function codes the name of a desired level and wavelet coefficient orientation into a string which is used by the 2D DWT functions to access and manipulate wavelet coefficients.
}
\details{
For the 1D wavelet transform (and others) the
\code{\link{accessC}} and \code{\link{accessD}} function extracts wavelet
coefficients from 1D wavelet decomposition objects.
For \code{\link{imwd.object}} class objects,
which are the 2D wavelet transforms of lattice objects (images)
the wavelet coefficients are stored within components of the list object
that underlies the imwd object.

This function provides an easy way to specify a resolution level and
orientation in a human readable way and this function then produces the
character string necessary to access the wavelet coefficients in an imwd object.

Note that this function \emph{does not} actually extract any coefficients itself.
}
\value{
A character string which codes the level and type of coefficients.
It reads wXLY X is the resolution level and Y is an integer corresponding to the orientation (1=horizontal, 2=vertical, 3=diagonal, 4=smoothed).
}
\seealso{
\code{\link{imwd}}, \code{\link{imwd.object}}
}
\examples{
#
# Generate the character string for the component of the imwd object
#
# The string associated with the diagonal detail at the third level...
#
lt.to.name(3, "DD")
# [1] "w3L3"
#
# Show how to access wavelet coefficients of imwd object.
#
# First, make up some data (using matrix/rnorm) and then subject it
# to an image wavelet transform.
#
tmpimwd <- imwd(matrix(rnorm(64),64,64))
#
# Get the horizontal coefficients at the 2nd level
#
tmpimwd[[ lt.to.name(2, "CD") ]]
# [1]  6.962251e-13  4.937486e-12  3.712157e-12 -3.064831e-12  6.962251e-13
# [6]  4.937486e-12  3.712157e-12 -3.064831e-12  6.962251e-13  4.937486e-12
# [11]  3.712157e-12 -3.064831e-12  6.962251e-13  4.937486e-12  3.712157e-12
# [16] -3.064831e-12
#
#
# If you want the coefficients returned as a matrix use the matrix function,
# i.e.
#
matrix(tmpimwd[[ lt.to.name(2, "CD") ]], 4,4)
#              [,1]          [,2]          [,3]          [,4]
#[1,]  6.962251e-13  6.962251e-13  6.962251e-13  6.962251e-13
#[2,]  4.937486e-12  4.937486e-12  4.937486e-12  4.937486e-12
#[3,]  3.712157e-12  3.712157e-12  3.712157e-12  3.712157e-12
#[4,] -3.064831e-12 -3.064831e-12 -3.064831e-12 -3.064831e-12
#
# Note that the dimensions of the matrix depend on the resolution level
# that you extract and dim = 2^level
}
\author{G P Nason}
\keyword{manip}
