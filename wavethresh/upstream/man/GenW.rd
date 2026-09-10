\name{GenW}
\alias{GenW}
\title{Generate (inverse) discrete wavelet transform matrix.}
\description{
This function generates a matrix that can perform the discrete wavelet transform (useful for understanding the DWT but use the fast algorithm coded in \code{\link{wd}} for general use). The function returns the matrix for the inverse transform. Since the matrix is orthogonal transpose the matrix to obtain the forward transform matrix. 
}
\usage{
GenW(n=8, filter.number=10, family="DaubLeAsymm", bc="periodic")
}
\arguments{
\item{n}{The order of the DWT matrix will be n times n. n should be a power of two.} 
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.} 
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".} 
\item{bc}{boundary conditions to use. This can be \code{periodic} or \code{symmetric} depending on whether you want the returned matrix to assume periodic or symmetric end-reflection boundary conditions.} 
}
\details{
The discrete wavelet transform is usually computed using the fast pyramid algorithm of Mallat. However, the transform can be written in a matrix form and this is useful for understanding what the fast transform does. One wouldn't normally use the matrix for performing the transform but use the fast transform function \code{\link{wd}} instead. 

The matrix returned by this function represents the inverse DWT. Since the matrix (and transform) is orthogonal one can obtain the matrix representation of the forward transform simply by transposing the matrix using the \code{t} function in S-Plus. 

The returned matrix is organised as follows. The first column always corresponds to the linear combination corresponding to the scaling function coefficient (so the column is constant. The next \code{n/2} columns correspond to the finest scale wavelet coefficients; the next \code{n/4} columns to the next finest scale and so on until the last column which corresponds to the coarsest scale wavelet coefficients. 

The matrix is computed by performing successive fast DWTs on unit vectors. 

}
\value{
A matrix of order \code{n} that contains the inverse discrete wavelet transform. 
}
\section{RELEASE}{Version 3.2 Copyright Guy Nason 1998 }
\seealso{
\code{\link{wd}}, \code{\link{wr}}. 
}
\examples{
#
# Generate the wavelet transform matrix corresponding to the Haar wavelet
# transform of order 8
#
haarmat <- GenW(8, filter.number=1, family="DaubExPhase")
#
# Let's look at this matrix
#
#haarmat
#          [,1]       [,2]       [,3]       [,4]       [,5] [,6] [,7]       [,8] 
#[1,] 0.3535534  0.7071068  0.0000000  0.0000000  0.0000000  0.5  0.0  0.3535534
#[2,] 0.3535534 -0.7071068  0.0000000  0.0000000  0.0000000  0.5  0.0  0.3535534
#[3,] 0.3535534  0.0000000  0.7071068  0.0000000  0.0000000 -0.5  0.0  0.3535534
#[4,] 0.3535534  0.0000000 -0.7071068  0.0000000  0.0000000 -0.5  0.0  0.3535534
#[5,] 0.3535534  0.0000000  0.0000000  0.7071068  0.0000000  0.0  0.5 -0.3535534
#[6,] 0.3535534  0.0000000  0.0000000 -0.7071068  0.0000000  0.0  0.5 -0.3535534
#[7,] 0.3535534  0.0000000  0.0000000  0.0000000  0.7071068  0.0 -0.5 -0.3535534
#[8,] 0.3535534  0.0000000  0.0000000  0.0000000 -0.7071068  0.0 -0.5 -0.3535534
#
# As noted above the first column is the l.c. corresponding to the scaling
# function coefficient and then the l.c.s corresponding to the wavelet
# coefficients from the finest to the coarsest.
#
# The above matrix represented the inverse DWT. Let's compute the forward
# transform matrix representation:
#
#t(haarmat)
#          [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]       [,8] 
#[1,] 0.3535534  0.3535534  0.3535534  0.3535534  0.3535534  0.3535534  0.3535534  0.3535534
#[2,] 0.7071068 -0.7071068  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#[3,] 0.0000000  0.0000000  0.7071068 -0.7071068  0.0000000  0.0000000  0.0000000  0.0000000
#[4,] 0.0000000  0.0000000  0.0000000  0.0000000  0.7071068 -0.7071068  0.0000000  0.0000000
#[5,] 0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.7071068 -0.7071068
#[6,] 0.5000000  0.5000000 -0.5000000 -0.5000000  0.0000000  0.0000000  0.0000000  0.0000000
#[7,] 0.0000000  0.0000000  0.0000000  0.0000000  0.5000000  0.5000000 -0.5000000 -0.5000000
#[8,] 0.3535534  0.3535534  0.3535534  0.3535534 -0.3535534 -0.3535534 -0.3535534 -0.3535534
#
#
}
\keyword{array}
\author{G P Nason}
