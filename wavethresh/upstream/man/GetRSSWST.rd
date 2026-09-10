\name{GetRSSWST}
\alias{GetRSSWST}
\title{Computes estimate of error for function estimate. }
\description{
Computes estimate of error for function estimate. Given noisy data and a threshold value this function uses Nason's 1996 two-fold cross-validation algorithm, but using packet ordered non-decimated wavelet transforms to compute two estimates of an underlying ``true'' function and uses them to compute an estimate of the error in estimating the truth. 
}
\usage{
GetRSSWST(ndata, threshold, levels, family = "DaubLeAsymm", 
	filter.number = 10, type = "soft", norm = l2norm, verbose = 0, 
	InverseType = "average")

}
\arguments{
\item{ndata}{the noisy data. This is a vector containing the signal plus noise. The length of this vector should be a power of two.} 
\item{threshold}{the value of the threshold that you wish to compute the error of the estimate at} 
\item{levels}{the levels over which you wish the threshold value to be computed (the threshold that is used in computing the estimate and error in the estimate). See the explanation for this argument in the \code{\link{threshold.wst}} function. }
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".} 
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.} 
\item{type}{whether to use hard or soft thresholding. See the explanation for this argument in the \code{\link{threshold.wst}} function.}
\item{norm}{which measure of distance to judge the dissimilarity between the estimates. The functions \code{\link{l2norm}} and \code{\link{linfnorm}} are suitable examples.} 
\item{verbose}{If \code{TRUE} then informative messages are printed during the progression of the function, otherwise they are not.} 
\item{InverseType}{The possible options are "average" or "minent". The former uses basis averaging to form estimates of the unknown function. The "minent" function selects a basis using the Coifman and Wickerhauser, 1992 algorithm to select a basis to invert.} 

}
\details{
This function implements the component of the cross-validation method detailed by Nason, 1996 for computing an estimate of the error between an estimate and the ``truth''. The difference here is that it uses the packet ordered non-decimated wavelet transform rather than the standard Mallat \code{\link{wd}} discrete wavelet transform. As such it is an examples of the translation-invariant denoising of Coifman and Donoho, 1995 but uses cross-validation to choose the threshold rather than SUREshrink. 

Note that the procedure outlined above can use \code{\link{AvBasis}} basis averaging or basis selection and inversion using the Coifman and Wickerhauser, 1992 best-basis algorithm 

}
\value{
A real number which is estimate of the error between estimate and truth at the given threshold. 
}
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1995 }
\seealso{
\code{\link{linfnorm}}, \code{\link{linfnorm}}, \code{\link{wstCV}}, \code{\link{wstCVl}}. 

}
\examples{
#
# This function performs the error estimation step for the
# \code{\link{wstCV}} function and so is not intended for
# user use. 
#


}
\keyword{manip}
\author{G P Nason}
