\name{first.last.dh}
\alias{first.last.dh}
\usage{
first.last.dh(LengthH, DataLength, type = "wavelet", bc = "periodic",
    firstk = c(0, DataLength - 1))
}
\arguments{
\item{LengthH}{The length of the smoothing (C) filter}
\item{DataLength}{The length of the data that you wish to transform}
\item{type}{The type of wavelet transform, \code{wavelet} or \code{station}
	for decimated and nondecimated transforms respectively.}
\item{bc}{Boundary conditions, \code{periodic} or \code{symmetric}}
\item{firstk}{The first k index, leave as default}

}
\title{Build special first/last database for some wavelet density functions}
\description{
This function builds a special first/last database for some of the
wavelet density estimation functions written by David Herrick and
described in his PhD thesis.

See \code{\link{first.last}} to see what this kind of function does.
}
\details{
Description says all.
}
\value{
A list with several components in exactly the same format as
for \code{\link{first.last}}.
}
\seealso{\code{\link{dencvwd}},\code{\link{first.last}},\code{\link{wd.dh}}}
\author{David Herrick}
\keyword{manip}
