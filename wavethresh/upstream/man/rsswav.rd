\name{rsswav}
\alias{rsswav}
\title{Compute mean residual sum of squares for odd prediction of even ordinates and vice versa}
\usage{
rsswav(noisy, value = 1, filter.number = 10, family = "DaubLeAsymm", 
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
Compute mean of residual sum of squares (RSS) for odd prediction of even
ordinates and vice versa using wavelet shrinkage with a specified threshold.
This is a subsidiary routine of the \code{\link{WaveletCV}}
cross validation function.
A version implemented in C exists called \code{\link{Crsswav}}.
}
\details{
\bold{Note}: a faster C based implementation of this function called
\code{\link{Crsswav}} is available.
It takes the same arguments and returns the same values.

Two-fold cross validation can be computed for a wd object using the "cv" policy option in \code{\link{threshold.wd}}.
As part of this procedure for each threshold value that the CV optimisation algorithm selects a RSS value must be computed (the CV optimisation algorithm seeks to minimize this RSS value).

The RSS value computed is this. First, the even and odd indexed values are separated. The even values are used to construct an estimate of the odd true values using wavelet shrinkage with the given threshold. The sum of squares between the estimate and the noisy odds is computed. An equivalent calculation is performed by swapping the odds and evens. The two RSS values are then averaged and the average returned. This algorithm is described more fully in Nason, (1996).
}
\value{
A list with the following components
\item{ssq}{The RSS value that was computed}
\item{df}{The dof value computed on the thresholded wavelet transform of the
data with the given threshold and thresholding options. (Although this is not really used for anything).}
\item{value}{The value argument that was specified.}
\item{type}{the \code{thresh.type} argument that was specified.}
\item{lev}{The vector \code{ll:(nlevelsWT(noisy)-1)} (i.e. the levels that were thresholded).}
}
\seealso{\code{\link{Crsswav}},\code{\link{threshold.wd}}, \code{\link{WaveletCV}}}
\author{G P Nason}
\keyword{smooth}
