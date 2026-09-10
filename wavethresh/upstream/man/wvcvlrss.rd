\name{wvcvlrss}
\alias{wvcvlrss}
\title{Computes estimate of error for function estimate. }
\description{
This function is merely a call to the \code{\link{GetRSSWST}} function. 
}
\usage{
wvcvlrss(threshold, ndata, levels, type, filter.number, family, norm, 
	verbose, InverseType)
}
\arguments{
\item{threshold}{the value of the threshold that you wish to compute the error of the estimate at}
\item{ndata}{the noisy data. This is a vector containing the signal plus noise. The length of this vector should be a power of two.}
\item{levels}{the levels over which you wish the threshold value to be computed (the threshold that is used in computing the estimate and error in the estimate). See the explanation for this argument in the \code{\link{threshold.wst}} function. }
\item{type}{whether to use hard or soft thresholding. See the explanation for this argument in the \code{\link{threshold.wst}} function.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{norm}{which measure of distance to judge the dissimilarity between the estimates. The functions \code{\link{l2norm}} and \code{\link{linfnorm}} are suitable examples.} 
\item{verbose}{If \code{TRUE} then informative messages are printed during the progression of the function, otherwise they are not.}
\item{InverseType}{The possible options are "average" or "minent". The former uses basis averaging to form estimates of the unknown function. The "minent" function selects a basis using the Coifman and Wickerhauser, 1992 algorithm to select a basis to invert.}
}
\details{
This function is merely a call to the \code{\link{GetRSSWST}} function with a few arguments interchanged. In particular, the first two arguments are interchanged. This is to make life easier for use with the \code{nlminb} function which expects the first argument of the function it is trying to optimise to be the variable that the function is optimised over. 
}
\value{
A real number which is estimate of the error between estimate and truth at the given threshold. 
}
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1995 }
\seealso{
\code{\link{GetRSSWST}}. 
}
\examples{
#
# This function performs the error estimation step for the
# wstCVl function and so is not intended for
# user use. 
#
}
\keyword{manip}
\author{G P Nason}
