\name{make.dwwt}
\alias{make.dwwt}
\title{Compute diagonal of the matrix WWT}
\description{
Computes the values which specify the covariance structure of complex-valued wavelet coefficients. 
}
\usage{
make.dwwt(nlevels, filter.number = 3.1, family = "LinaMayrand") 
}
\arguments{
\item{nlevels}{The number of levels of the wavelet decomposition.}
\item{filter.number, family}{Specifies the wavelet used; see filter.select for more details.} 
}
\details{
If real-valued signals are decomposed by a discrete wavelet transform using a complex-valued Daubechies wavelet (as described by Lina & Mayrand (1995)), the resulting coefficients are complex-valued. The covariance structure of these coefficients are determined by the diagonal entries of the matrix
\eqn{WW^T}. This function computes these values for use in shrinkage. For more details, see Barber & Nason (2004) 
}
\value{
A vector giving the diagonal elements of \eqn{WW^T}. 
}
\section{RELEASE}{
Part of the CThresh addon to WaveThresh. Copyright Stuart Barber and Guy Nason 2004. }
\seealso{
\code{\link{cthresh}}
}
\keyword{manip}
\author{Stuart Barber}
