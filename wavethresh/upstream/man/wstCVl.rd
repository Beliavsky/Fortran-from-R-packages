\name{wstCVl}
\alias{wstCVl}
\title{Performs two-fold cross-validation estimation using packet-ordered non-decimated wavelet transforms and a (vector) level-dependent threshold. }
\description{
Performs Nason's 1996 two-fold cross-validation estimation using packet-ordered non-decimated wavelet transforms and a (vector) level-dependent threshold.
}
\usage{
wstCVl(ndata, ll = 3, type = "soft", filter.number = 10, family = "DaubLeAsymm",
	tol = 0.01, verbose = 0, plot.it = FALSE, norm = l2norm, InverseType = "average",
	uvdev = madmad)
}
\arguments{
\item{ndata}{the noisy data. This is a vector containing the signal plus noise. The length of this vector should be a power of two.}
\item{ll}{the primary resolution for this estimation. Note that the primary resolution is \emph{problem-specific}: you have to find out which is the best value.}
\item{type}{whether to use hard or soft thresholding. See the explanation for this argument in the \code{\link{threshold.wst}} function.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{tol}{the cross-validation tolerance which decides when an estimate is sufficiently close to the truth (or estimated to be so).}
\item{verbose}{If \code{TRUE} then informative messages are printed during the progression of the function, otherwise they are not.}
\item{plot.it}{Whether or not to produce a plot indicating progress.}
\item{norm}{which measure of distance to judge the dissimilarity between the estimates. The functions \code{\link{l2norm}} and \code{\link{linfnorm}} are suitable examples.}
\item{InverseType}{The possible options are "average" or "minent". The former uses basis averaging to form estimates of the unknown function. The "minent" function selects a basis using the Coifman and Wickerhauser, 1992 algorithm to select a basis to invert.}
\item{uvdev}{Universal thresholding is used to generate an upper bound for the ideal threshold. This argument provides the function that computes an estimate of the variance of the noise for use with the universal threshold calculation (see \code{\link{threshold.wst}}).}
}
\details{
This function implements a modified version of the cross-validation method detailed by Nason, 1996 for computing an estimate of the error between an estimate and the ``truth''. The difference here is that it uses the packet ordered non-decimated wavelet transform rather than the standard Mallat wd discrete wavelet transform. As such it is an examples of the translation-invariant denoising of Coifman and Donoho, 1995 but uses cross-validation to choose the threshold rather than SUREshrink. 

Further, this function computes level-dependent thresholds. That is, it can compute a different threshold for each resolution level. 

Note that the procedure outlined above can use \code{\link{AvBasis}} basis averaging or basis selection and inversion using the Coifman and Wickerhauser, 1992 best-basis algorithm 
}
\value{
A list returning the results of the cross-validation algorithm. The list includes the following components: 

\item{ndata}{a copy of the input noisy data}
\item{xvwr}{a reconstruction of the best estimate computed using this algorithm. It is the inverse (computed depending on what the InverseType argument was) of the \code{xvwrWSTt} component.} 
\item{xvwrWSTt}{a thresholded version of the packet-ordered non-decimated wavelet transform of the noisy data using the best threshold discovered by this cross-validation algorithm.} 
\item{uvt}{the universal threshold used as the upper bound for the algorithm that tries to discover the optimal cross-validation threshold. The lower bound is always zero.}
\item{xvthresh}{the best threshold as discovered by cross-validation. Note that this is vector, a level-dependent threshold with one threshold value for each resolution level. The first entry corresponds to level \code{ll}, the last entry corresponds to level \code{nlevelsWT(ndata)-1} and the entries in between linearly to the levels in between. The \code{\link{wstCV}} function should be used to compute a global threshold.}
\item{optres}{The results from performing the optimisation using the \code{nlminb} function from Splus. This object contains many interesting components with information about how the optimisation went. See the \code{nlminb} help page for information.}
} 
\section{RELEASE}{
Version 3.6 Copyright Guy Nason 1995 
}
\seealso{
\code{\link{GetRSSWST}}, \code{\link{linfnorm}}, \code{\link{linfnorm}}, \code{\link{threshold.wst}}, \code{\link{wst}}, \code{\link{wst.object}}, \code{\link{wstCV}}
}
\examples{
#
# Example PENDING
#
}
\keyword{math}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
