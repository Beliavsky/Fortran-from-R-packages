\name{TrenchInverse}
\alias{TrenchInverse}
\title{compute the matrix inverse of a positive-definite Toepliz matrix }

\description{
The Trench algorithm (Golub and Vanload, 1983) is implemented in C and
interfaced to R.
This provides an expedient method for obtaining the
matrix inverse of the covariance matrix of n successive observations
from a stationary time series.
Some applications of this are discussed by McLeod and Krougly (2005).
}

\usage{TrenchInverse(G)}

\arguments{
  \item{G}{ a positive definite Toeplitz matrix }
}

\value{the matrix inverse of G is computed}

\references{ 
Golub, G. and Van Loan (1983).
Matrix Computations, 2nd Ed.
John Hoptkins University Press, Baltimore.
Algorithm 5.7-3.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.

}
\author{ A.I. McLeod }

\note{
TrenchInverse(x) assumes that x is a symmetric Toeplitz matrix but it
does not specifically test for this.
Instead it merely takes the first row of x and passes this directly
to the C code program which uses this more compact storage format.
The C code program then computes the inverse. 
An error message is given if the C code algorithm encounters
a non-positive definite input.
}

\section{Warning }{
You should test the input x using is.toeplitz(x) if you are not sure if
x is a symmetric Toeplitz matix.
}

\seealso{
  \code{\link{TrenchLoglikelihood}}, \code{\link{is.toeplitz}}, 
  \code{\link{DLLoglikelihood}}, 
  \code{\link{TrenchMean}}, \code{\link{solve}}
}

\examples{
#compute inverse of matrix and compare with result from solve
data(LakeHuron)
r<-acf(LakeHuron, plot=FALSE, lag.max=4)$acf
R<-toeplitz(c(r))
Ri<-TrenchInverse(R)
Ri2<-solve(R)
Ri
Ri2

#invert a matrix of order n and compute the maximum absolute error
# in the product of this inverse with the original matrix
n<-5	   
r<-0.8^(0:(n-1))
G<-toeplitz(r)
Gi<-TrenchInverse(G)
GGi<-crossprod(t(G),Gi)
id<-matrix(0, nrow=n, ncol=n)
diag(id)<-1
err<-max(abs(id-GGi))
err
}
\keyword{ts }
\keyword{array }

