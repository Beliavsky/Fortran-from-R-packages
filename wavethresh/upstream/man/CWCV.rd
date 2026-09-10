\name{CWCV}
\alias{CWCV}
\title{C Wavelet Cross-validation}
\description{
Two-fold wavelet shrinkage cross-validation (in C) 
}
\usage{
CWCV(ynoise, ll, x = 1:length(ynoise), filter.number = 10, family = 
	"DaubLeAsymm", thresh.type = "soft", tol = 0.01,
	maxits=500, verbose = 0, 
	plot.it = TRUE, interptype = "noise")

}
\arguments{
\item{ynoise}{A vector of dyadic (power of two) length that contains the noisy data that you wish to apply wavelet shrinkage by cross-validation to.}
\item{ll}{The primary resolution that you wish to assume. No wavelet coefficients that are on coarser scales than ll will be thresholded.}
\item{x}{This function is capable of producing informative plots. It can be useful to supply the x values corresponding to the ynoise values. Further this argument is returned by this function which can be useful for later processors.}
\item{filter.number}{This selects the smoothness of wavelet that you want to perform wavelet shrinkage by cross-validation.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{thresh.type }{this option specifies the thresholding type which can be "hard" or "soft".}
\item{tol}{this specifies the convergence tolerance for the cross-validation optimization routine (a golden section search).}
\item{maxits}{maximum number of iterations for the cross-validation optimization routine (a golden section search).}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default}
\item{plot.it}{If this is TRUE then plots of the universal threshold (used to obtain an upper bound on the cross-validation threshold) reconstruction and the resulting cross-validation estimate are produced.}
\item{interptype}{Can take two values noise or normal. This option controls how cross-validation compares the estimate formed by leaving out the data with the "left-out" data. If interptype="noise" then two noisy values are averaged to compare with the estimated curve in between, otherwise if interptype="normal" then the curve estimate is averaged either side of a noisy left-out point.}
}
\details{
Compute the two-fold cross-validated wavelet shrunk estimate given the noisy data ynoise according to the description given in Nason, 1996. 

You must specify a primary resolution given by \code{ll}. This must be specified individually on each data set and can itself be estimated using cross-validation (although I haven't written the code to do this). 

\bold{Note}. The two-fold cross-validation method performs very badly if the input data is correlated. In this case I would advise using the methods proposed in Donoho and Johnstone, 1995 or Johnstone and Silverman, 1997 which can be carried out in WaveThresh using the \code{\link{threshold}} function using the \code{policy="sure"} option. 

}
\value{
A list with the following components 

\item{x}{This is just the x that was input. It gets passed through more or less for convenience for the user.}
\item{ynoise}{A copy of the input ynoise noisy data.}
\item{xvwr}{The cross-validated wavelet shrunk estimate.}
\item{yuvtwr}{The universal thresholded version (note this is merely a starting point for the cross-validation algorithm. It should not be ta    ken seriously as an estimate. In particular its estimate of variance is likely to be inflated.) }
\item{xvthresh}{The cross-validated threshold}
\item{xvdof}{The number of non-zero coefficients in the cross-validated shrunk wavelet object (which is not returned).}
\item{uvdof}{The number of non-zero coefficients in the universal threshold shrunk wavelet object (which also is not returned)}
\item{xkeep}{always returns NULL!}
\item{fkeep}{always returns NULL!}

}
\note{Plots of the universal and cross-validated shrunk estimates might be plotted if \code{plot.it=TRUE.}}
\section{RELEASE}{Version 3.0 Copyright Guy Nason 1994 }
\seealso{
\code{\link{threshold}}. \code{\link{threshold.wd}}. 
}
\examples{
#
# This function is best used via the policy="cv" option in
# the threshold.wd function.
# See examples there.
#
}
\keyword{smooth}
\author{G P Nason}
