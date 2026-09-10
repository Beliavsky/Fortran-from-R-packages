\name{Chires6}
\alias{Chires6}
\title{Subsid routine for denproj (calcs scaling function coefs with cov)}
\usage{
Chires6(x, tau=1, J, filter.number=10, family="DaubLeAsymm", nT=20)
}
\arguments{
\item{x}{The data (random sample for density estimation)}
\item{tau}{Fine tuning parameter}
\item{J}{Resolution level}
\item{filter.number}{The smoothness of the wavelet, see \code{\link{filter.select}}}
\item{family}{The family of the wavelet, see \code{\link{family}}}
\item{nT}{The number of iterations in the Daubechies-Lagarias algorithm}
}
\description{Function is essentially the same as \code{\link{Chires5}}
but also returns covariances between coefficients.
A subsidiary routine for \code{\link{denproj}}. Not intended for direct
user use. 
}
\details{
As description
}
\value{
A list with the following components:
\item{coef}{The scaling function coefficients}
\item{covar}{The coefficients' covariance matrix}
\item{klim}{The integer translates of the scaling functions used}
\item{p}{The primary resolution, calculated in code as tau*2^J}
\item{filter}{The usual filter information, see \code{\link{filter.select}}}
\item{n}{The length of the data \code{x}}
\item{res}{A list containing components: \code{p}, as above, \code{tau}
	as input and \code{J} as above. This summarizes the resolution
	information}
}
\seealso{\code{\link{Chires6}},\code{\link{denproj}}}
\author{David Herrick}
\keyword{smooth}
