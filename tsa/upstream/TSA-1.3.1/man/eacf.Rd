\name{eacf}
\alias{eacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the sample extended acf (ESACF)}
\description{
Computes the sample extended acf (ESACF) for the time series stored in z.
The matrix of ESACF with the AR order up to ar.max and the MA order
up to ma.max is stored in the matrix EACFM.
}
\usage{
eacf(z, ar.max = 7, ma.max = 13)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{the time series data}
  \item{ar.max}{maximum AR order; default=7}
  \item{ma.max}{maximum MA order; default=13}
}
\value{
A list containing the following two components:  
\item{eacf}{a matrix of sample extended ACF} 
\item{symbol}{corresponding matrix of symbols indicating the significance of
the ESACF} 
Side effect of the eacf function:
The function prints a coded ESACF table with
significant values denoted by * and nosignificant values by 0. 
}
\references{Tsay, R. and Tiao, G. (1984). "Consistent Estimates of Autoregressive Parameters and
Extended Sample Autocorrelation Function for Stationary and Nonstationary
ARMA Models." Journal of the American Statistical Association, 79 (385), pp.
84-96.}
\author{Kung-Sik Chan}
\examples{
data(arma11.s)
eacf(arma11.s)
}
\keyword{methods}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
