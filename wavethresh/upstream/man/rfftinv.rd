\name{rfftinv}
\alias{rfftinv}
\title{Inverse real FFT, inverse of rfft}
\usage{
rfftinv(rz, n = length(rz))
}
\arguments{
\item{rz}{The Fourier coefficients to invert}
\item{n}{The number of coefficients}
}
\description{
Inverse function of \code{\link{rfft}}
}
\details{
Just the inverse function of \code{\link{rfft}}.
}
\value{
The inverse FT of the input
}
\seealso{\code{\link{rfft}}}
\author{Bernard Silverman}
\keyword{math}
