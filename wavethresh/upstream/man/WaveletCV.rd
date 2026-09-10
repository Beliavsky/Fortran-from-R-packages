\name{WaveletCV}
\alias{WaveletCV}
\title{Wavelet cross-validation}
\usage{
WaveletCV(ynoise, x = 1:length(ynoise), filter.number = 10, family = 
	"DaubLeAsymm", thresh.type = "soft", tol = 0.01, verbose = 0, 
	plot.it = TRUE, ll=3)
}
\arguments{
\item{ynoise}{A vector of dyadic (power of two) length that contains the
	noisy data that you wish to apply wavelet shrinkage by
	cross-validation to.}
\item{x}{This function is capable of producing informative plots.
	It can be useful to supply the x values corresponding to the
	\code{ynoise} values.
	Further this argument is returned by this function which can
	be useful for later processors.}
\item{filter.number}{This selects the smoothness of wavelet that you want to perform wavelet shrinkage by cross-validation.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{thresh.type}{this option specifies the thresholding type which can be "hard" or "soft".}
\item{tol}{this specifies the convergence tolerance for the cross-validation optimization routine (a golden section search).}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
\item{plot.it}{If this is TRUE then plots of the universal threshold (used to obtain an upper bound on the cross-validation threshold) reconstruction and the resulting cross-validation estimate are produced.}
\item{ll}{The primary resolution that you wish to assume. No wavelet coefficients that are on coarser scales than ll will be thresholded.}
}
\description{
Two-fold wavelet shrinkage cross-validation (there is a faster C based
version \code{\link{CWCV}}.)
}
\details{
\bold{Note:} a faster C based implementation of this function called
\code{\link{CWCV}} is available.
It takes the same arguments (although it has one extra minor argument) and returns the same values.

Compute the two-fold cross-validated wavelet shrunk estimate given the noisy data ynoise according to the description given in Nason, 1996.

You must specify a primary resolution given by \code{ll}. This must be specified individually on each data set and can itself be estimated using cross-validation (although I haven't written the code to do this).

Note. The two-fold cross-validation method performs very badly if the input data is correlated. In this case I would advise using other methods.
}
\value{
A list with the following components
\item{x}{This is just the x that was input. It gets passed through more or less for convenience for the user.}
\item{ynoise}{A copy of the input ynoise noisy data.}
\item{xvwr}{The cross-validated wavelet shrunk estimate.}
\item{yuvtwr}{The universal thresholded version (note this is merely a starting point for the cross-validation algorithm. It should not be taken seriously as an estimate. In particular its estimate of variance is likely to be inflated.)}
\item{xvthresh}{The cross-validated threshold}
\item{uvthresh}{The universal threshold (again, don't take this value too seriously. You might get better performance using the threshold function directly with specialist options.}
\item{xvdof}{The number of non-zero coefficients in the cross-validated shrunk wavelet object (which is not returned).}
\item{uvdof}{The number of non-zero coefficients in the universal threshold shrunk wavelet object (which also is not returned)}
\item{xkeep}{always returns NULL!}
\item{fkeep}{always returns NULL!}
}
\seealso{\code{\link{CWCV}},\code{\link{Crsswav}},\code{\link{rsswav}},\code{\link{threshold.wd}}}
\examples{
#
# This function is best used via the policy="cv" option in
# the threshold.wd function.
# See examples there.
#
}
\author{G P Nason}
\keyword{smooth}
