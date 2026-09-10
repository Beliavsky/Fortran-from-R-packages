\name{ipndacw}
\alias{ipndacw}
\title{Compute inner product matrix of discrete non-decimated autocorrelation wavelets.}
\description{
This function computes the inner product matrix of discrete non-decimated autocorrelation wavelets. 
}
\usage{
ipndacw(J, filter.number = 10, family = "DaubLeAsymm", tol = 1e-100, verbose
         = FALSE, \dots)
}
\arguments{
\item{J}{Dimension of inner product matrix required. This number should be a negative integer.}
\item{filter.number}{The index of the wavelet used to compute the inner product matrix.} 
\item{family}{The family of wavelet used to compute the inner product matrix.} 
\item{tol}{In the brute force computation for Daubechies compactly supported wavelets many inner product computations are performed. This tolerance discounts any results which are smaller than \code{tol} which effectively defines how long the inner product/autocorrelation products are.}
\item{verbose}{If \code{TRUE} then informative messages are printed. Some of these can be quite fun as the function tells you whether precomputed matrices are being used, how much computation needs to be done and so forth. }
\item{\dots}{any other arguments}
}
\details{
This function computes the inner product matrix of the discrete non-decimated autocorrelation wavelets. This matrix is used to correct the wavelet periodogram as a step to turning it into a evolutionary wavelet spectral estimate. The matrix returned by ipndacw is the one called A in the paper by Nason, von Sachs and Kroisandt. 

For the Haar wavelet the matrix is computed by using the analytical formulae in the paper by Nason, von Sachs and Kroisandt and is hence very fast and efficient and can be used for large values of -J. 

For other Daubechies compactly supported wavelets the matrix is computed directly by autocorrelating discrete non-decimated wavelets at different scales and then forming the inner products of these. A function that computes the autocorrelation wavelets themselves is \code{\link{PsiJ}}. This \emph{brute force} computation is slow and memory inefficient hence \code{ipndacw} contains a mechanism that stores any inner product matrix that it creates according to a naming scheme defined by the convention defined in \code{\link{rmname}}. The stored matrices are assigned to the user-visible environment \code{\link{WTEnv}}.

These stored matrices can be used in future computations by the following automatic procedure: 
\describe{
\item{1}{The \code{\link{rmget}} looks to see whether previous computations have been performed that might be useful.}
\item{2}{If a matrix of higher order is discovered then the appropriate top-left submatrix is returned, otherwise...}
\item{3}{If the right order of matrix is found it is returned, otherwise ...}
\item{4}{If a matrix of \emph{smaller} order is found it is used as the top-left submatrix of the answer. The remaining elements to the right of and below the submatrix are computed and then the whole matrix is returned, otherwise...} 
\item{5}{If none are found then the whole matrix is computed in C and returned.}
}
In this way a particular matrix for a given wavelet need only be computed once. 

}
\value{
A matrix of order (-J)x(-J) containing the inner product matrix of the discrete non-decimated autocorrelation matrices. 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern. }
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{ewspec}}, \code{\link{PsiJ}}, \code{\link{rmname}}, \code{\link{rmget}}, \code{\link{filter.select}}. 
}
\examples{
#
# Let us create the 4x4 inner product matrix for the Haar wavelet.
# We'll turn on the jolly verbose messages as well. 
#
ipndacw(-4, filter.number=1, family="DaubExPhase", verbose=TRUE)
#Computing ipndacw
#Calling haarmat
#Took  0.0699999  seconds
#       -1     -2     -3     -4 
#-1 1.5000 0.7500 0.3750 0.1875
#-2 0.7500 1.7500 1.1250 0.5625
#-3 0.3750 1.1250 2.8750 2.0625
#-4 0.1875 0.5625 2.0625 5.4375
#
# If we do this again it will use the precomputed version
#
ipndacw(-4, filter.number=1, family="DaubExPhase", verbose=TRUE)
#Computing ipndacw
#Returning precomputed version: using  4 
#Took  0.08  seconds
#       -1     -2     -3     -4 
#-1 1.5000 0.7500 0.3750 0.1875
#-2 0.7500 1.7500 1.1250 0.5625
#-3 0.3750 1.1250 2.8750 2.0625
#-4 0.1875 0.5625 2.0625 5.4375
#
# Let's use a smoother wavelet from the least-asymmetric family
# and generate the 6x6 version.
#
ipndacw(-6, filter.number=10, family="DaubLeAsymm", verbose=TRUE)
#Computing ipndacw
#Took  0.95  seconds
#             -1           -2           -3           -4           -5 
#-1 1.839101e+00 3.215934e-01 4.058155e-04 8.460063e-06 4.522125e-08
#-2 3.215934e-01 3.035353e+00 6.425188e-01 7.947454e-04 1.683209e-05
#-3 4.058155e-04 6.425188e-01 6.070419e+00 1.285038e+00 1.589486e-03
#-4 8.460063e-06 7.947454e-04 1.285038e+00 1.214084e+01 2.570075e+00
#-5 4.522125e-08 1.683209e-05 1.589486e-03 2.570075e+00 2.428168e+01
#-6 5.161675e-10 8.941666e-08 3.366416e-05 3.178972e-03 5.140150e+00
#             -6 
#-1 5.161675e-10
#-2 8.941666e-08
#-3 3.366416e-05
#-4 3.178972e-03
#-5 5.140150e+00
#-6 4.856335e+01
#
}
\keyword{algebra}
\author{G P Nason}
