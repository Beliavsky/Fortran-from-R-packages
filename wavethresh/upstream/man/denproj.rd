\name{denproj}
\alias{denproj}
\title{Calculate empirical scaling function coefficients of a p.d.f.}
\usage{
denproj(x, tau=1, J, filter.number=10, family="DaubLeAsymm", covar=FALSE, nT=20)
}
\arguments{
\item{x}{Vector containing the data. This can be of any length.}
\item{J}{The resolution level at which the empirical scaling function
	coefficients are to be calculated.}
\item{tau}{This parameter allows non-dyadic resolutions to be used,
	since the resolution is specified as \code{tau * 2J}.}
\item{filter.number}{The filter number of the wavelet basis to be used.}
\item{family}{The family of wavelets to use, can be "DaubExPhase" or
"DaubLeAsymm".}
\item{covar}{Logical variable. If TRUE then covariances of the empirical
	scaling function coefficients are also calculated.}
\item{nT}{The number of iterations to be performed in the
Daubechies-Lagarias algorithm,
which is used to evaluate the scaling functions of the specified wavelet basis at the data points.}
}
\description{
Calculates empirical scaling function coefficients of the probability
density function from sample of data from that density,
usually at some "high" resoloution.
}
\details{
This projection of data onto a high resolution wavelet space is described in
detail in Chapter 3 of Herrick (2000).
The maximum and minimum values of \code{k} for which the empirical scaling
function coefficient is non-zero are determined and
the coefficients calculated for all k between these limits as
\code{sum(phiJk(xi))/n}.
The scaling functions are evaluated at the data points efficiently,
using the Daubechies-Lagarias algorithm (Daubechies & Lagarias (1992)).
Coded kindly by Brani Vidakovic.

Herrick, D.R.M. (2000) Wavelet Methods for Curve and Surface Estimation. PhD Thesis, University of Bristol.

Daubechies, I. & Lagarias, J.C. (1992). Two-Scale Difference Equations II. Local Regularity, Infinite Products of Matrices and Fractals. SIAM Journal on Mathematical Analysis, 24(4), 1031--1079.
}
\value{
A list with components:
\item{coef}{A vector containing the empirical scaling function coefficients. This starts with the first non-zero coefficient, ends with the last non-zero coefficient and contains all coefficients, including zeros, in between.}
\item{covar}{Matrix containing the covariances, if requested.}
\item{klim}{The maximum and minimum values of k for which the empirical
	scaling function coefficients cJk are non-zero.}
\item{p}{The primary resolution \code{tau * 2J}.}
\item{filter}{A list containing the filter.number and family specified inthe function call.}
\item{n}{The length of the data vector x.}
\item{res}{A list containing the values of \code{p}, \code{tau} and \code{J}.}
}
\author{David Herrick}
\seealso{\code{\link{Chires5}}, \code{\link{Chires6}}, \code{\link{denwd}},
	\code{\link{denwr}}}
\examples{
# Simulate data from the claw density and find the 
# empirical scaling function coefficients
data <- rclaw(100)
datahr <- denproj(data, J=8, filter.number=4,family="DaubLeAsymm")
}
\keyword{smooth}
