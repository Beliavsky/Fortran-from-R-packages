\name{Crsswav}
\alias{Crsswav}
\title{Wrapper to C code version of rsswav}
\usage{
Crsswav(noisy, value = 1, filter.number = 10, family = "DaubLeAsymm",
    thresh.type = "hard", ll = 3)
}
\arguments{
\item{noisy}{A vector of dyadic (power of two) length that contains the noisy data that you wish to compute the averaged RSS for.}
\item{value}{The specified threshold.}
\item{filter.number}{This selects the smoothness of wavelet that you want to perform wavelet shrinkage by cross-validation.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{thresh.type}{this option specifies the thresholding type which can be "hard" or "soft".}
\item{ll}{The primary resolution that you wish to assume. No wavelet coefficients that are on coarser scales than ll will be thresholded.}
}
\description{
Crsswav is called by \code{\link{WaveletCV}} which is itself
called by \code{\link{threshold.wd}} to carry out its cross-validation
policy.
}
\details{
Description says all
}
\value{
Same value as for \code{\link{rsswav}}
}
\seealso{\code{\link{rsswav}}, \code{\link{WaveletCV}}}
\author{G P Nason}
\keyword{smooth}
