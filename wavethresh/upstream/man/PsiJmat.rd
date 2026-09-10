\name{PsiJmat}
\alias{PsiJmat}
\title{Compute discrete autocorrelation wavelets but return result in matrix form.}
\description{
This function computes discrete autocorrelation wavelets using the \code{\link{PsiJ}} function but it returns the results as a matrix rather than a list object.
}
\usage{
PsiJmat(J, filter.number = 10, family = "DaubLeAsymm", OPLENGTH=10^7)
}
\arguments{
\item{J}{Discrete autocorrelation wavelets will be computed for scales -1 up to scale J. This number should be a negative integer.}
\item{filter.number}{The index of the wavelet used to compute the discrete autocorrelation wavelets.}
\item{family}{The family of wavelet used to compute the discrete autocorrelation wavelets.} 
\item{OPLENGTH}{This integer variable defines some workspace of length OPLENGTH. The code uses this workspace. If the workspace is not long enough then the routine will stop and probably tell you what OPLENGTH should be set to.}
}
\details{
The discrete autocorrelation wavelet values are computed using the \code{\link{PsiJ}} function. This function merely organises them into a matrix form. 
}
\value{
A matrix containing -J rows and a number of columns less than OPLENGTH. Each row contains the values of the discrete autocorrelation wavelet for a different scale. Row one contains the scale -1 coefficients, row two contains the scale -2, and so on. 

The number of columns is an odd number. The middle position of each row is the value of the discrete autocorrelation wavelet at zero --- this is always 1. The discrete autocorrelation wavelet is symmetric about this point. 

\emph{Important} Apart from the central element none of the other columns line up in this way. This could be improved upon. 

}
\references{
Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern. }
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{PsiJ}}
}
\examples{
#
# As a simple first examples we shall compute the matrix containing
# the discrete autocorrelation wavelets up to scale 3.
#
PsiJmat(-3, filter.number=1, family="DaubExPhase")
#Computing PsiJ
#Took  0.25  seconds
#       [,1]  [,2]   [,3] [,4]   [,5]  [,6]   [,7] [,8]   [,9] [,10]  [,11] 
#[1,]  0.000  0.00  0.000  0.0  0.000  0.00 -0.500    1 -0.500  0.00  0.000
#[2,]  0.000  0.00  0.000  0.0 -0.250 -0.50  0.250    1  0.250 -0.50 -0.250
#[3,] -0.125 -0.25 -0.375 -0.5 -0.125  0.25  0.625    1  0.625  0.25 -0.125
#     [,12]  [,13] [,14]  [,15] 
#[1,]   0.0  0.000  0.00  0.000
#[2,]   0.0  0.000  0.00  0.000
#[3,]  -0.5 -0.375 -0.25 -0.125
#
# Note that this contains 3 rows (since J=-3).
# Each row contains the same discrete autocorrelation wavelet at different
# scales and hence different resolutions.
# Compare to the output given by PsiJ for the
# equivalent wavelet and scales.
# Note also that apart from column 8 which contains 1 (the value of the
# ac wavelet at zero) none of the other columns line up. E.g. the value of
# this wavelet at 1/2 is -0.5: this appears in columns 9, 10 and 12
# we could have written it differently so that they should line up.
# I might do this in the future. 
#
#
# Let's compute the matrix containing the discrete autocorrelation
# wavelets up to scale 6 using Daubechies N=10 least-asymmetric
# wavelets.
#
P6mat <- PsiJmat(-6, filter.number=10, family="DaubLeAsymm")
#
# What is the dimension of this matrix?
#
dim(P6mat)
#[1]    6 2395
#
# Hmmm. Pretty large, so we shan't print it out.
#
# However, these are the ac wavelets... Therefore if we compute their
# inner product we should get the same as if we used the ipndacw
# function directly.
#
P6mat %*% t(P6mat)
#             [,1]         [,2]         [,3]         [,4]         [,5] 
#[1,] 1.839101e+00 3.215934e-01 4.058155e-04 8.460063e-06 4.522125e-08
#[2,] 3.215934e-01 3.035353e+00 6.425188e-01 7.947454e-04 1.683209e-05
#[3,] 4.058155e-04 6.425188e-01 6.070419e+00 1.285038e+00 1.589486e-03
#[4,] 8.460063e-06 7.947454e-04 1.285038e+00 1.214084e+01 2.570075e+00
#[5,] 4.522125e-08 1.683209e-05 1.589486e-03 2.570075e+00 2.428168e+01
#[6,] 5.161675e-10 8.941666e-08 3.366416e-05 3.178972e-03 5.140150e+00
#             [,6] 
#[1,] 5.161675e-10
#[2,] 8.941666e-08
#[3,] 3.366416e-05
#[4,] 3.178972e-03
#[5,] 5.140150e+00
#[6,] 4.856335e+01
#
# Let's check it against the ipndacw call
#
ipndacw(-6, filter.number=10, family="DaubLeAsymm")
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
# Yep, they're the same.
#
}
\keyword{manip}
\author{G P Nason}
