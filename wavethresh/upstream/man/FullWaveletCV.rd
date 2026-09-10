\name{FullWaveletCV}
\alias{FullWaveletCV}
\title{Perform whole wavelet cross-validation in C code}
\usage{
FullWaveletCV(noisy, ll = 3, type = "soft", filter.number = 10, family =
    "DaubLeAsymm", tol = 0.01, verbose = 0)
}
\arguments{
\item{noisy}{A vector of dyadic (power of two) length that contains the
	noisy data that you wish to apply wavelet shrinkage by
	cross-validation to.}
\item{ll}{The primary resolution that you wish to assume. No wavelet coefficients that are on coarser scales than ll will be thresholded.}
\item{type}{this option specifies the thresholding type which can be "hard" or "soft".}
\item{filter.number}{This selects the smoothness of wavelet that you want to perform wavelet shrinkage by cross-validation.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{tol}{this specifies the convergence tolerance for the cross-validation optimization routine (a golden section search).}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
}
\description{
Perform whole wavelet cross-validation in C code. This
routine equivalent to \code{\link{CWCV}} except that
more preparatory material is passed to C code for speed. 

The major difference is that \bold{only} the cross-validated wavelet
threshold is returned.
}
\details{
Description says all
}
\value{
The cross-validated wavelet threshold.
}
\seealso{\code{\link{CWCV}}}
\author{G P Nason}
\keyword{smooth}
