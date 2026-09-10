\name{TOgetthrda1}
\alias{TOgetthrda1}
\alias{TOgetthrda2}
\alias{TOkolsmi.chi2}
\alias{TOonebyone1}
\alias{TOonebyone2}
\alias{TOshrinkit}
\title{Subsidiary routines for Ogden and Parzen's wavelet shrinkage methods}
\usage{
TOgetthrda1(dat, alpha)
TOgetthrda2(dat, alpha)
TOkolsmi.chi2(dat)
TOonebyone1(dat, alpha)
TOonebyone2(dat, alpha)
TOshrinkit(coeffs, thresh)
}
\arguments{
\item{dat}{data}
\item{alpha}{a p-value, generally smoothing parameter}
\item{coeffs}{Some coefficients to be shrunk}
\item{thresh}{a threshold}
}
\description{
Corresponds to the wavelet thresholding routine developed by
Ogden and Parzen (1994) Data dependent wavelet thresholding in nonparametric
regression with change-point applications. \emph{Tech Rep 176},
University of South Carolina, Department of Statistics.
}
\details{
Not intended for direct use.
}
\value{
Various depending on the function
}
\seealso{\code{\link{TOthreshda1}},\code{\link{TOthreshda2}},\code{\link{threshold}}}
\author{Todd Ogden}
\keyword{smooth}
