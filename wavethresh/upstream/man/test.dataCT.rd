\name{test.dataCT}
\alias{test.dataCT}
\title{Test functions for wavelet regression and thresholding }
\description{
This function evaluates the "blocks", "bumps", "heavisine" and "doppler" test functions of Donoho & Johnstone (1994b) and the piecewise polynomial test function of Nason & Silverman (1994). The function also generates data sets consisting of the specified function plus uncorrelated normally distributed errors. 
}
\usage{
test.dataCT(type = "ppoly", n = 512, signal = 1, rsnr = 7, plotfn = FALSE)
}
\arguments{
\item{type}{Test function to be computed. Available types are "ppoly" (piecewise polynomial), "blocks", "bumps", "heavi" (heavisine), and "doppler".}
\item{n}{Number of equally spaced data points on which the function is evaluated. }
\item{signal}{Scaling parameter; the function will be scaled so that the standard deviation of the data points takes this value.}
\item{rsnr}{Root signal-to-noise ratio. Specifies the ratio of the standard deviation of the function to the standard deviation of the simulated errors.}
\item{plotfn}{If \code{plotfn=TRUE}, then the test function and the simulated data set are plotted} 
}
\value{
A list with the following components: 
\item{x}{The points at which the test function is evaluated.}
\item{y}{The values taken by the test function.}
\item{ynoise}{The simulated data set.}
\item{type}{The type of function generated, identical to the input parameter type.}
\item{rsnr}{The root signal-to-noise ratio of the simulated data set, identical to the input parameter rsnr.} 
}
\section{Side effects}{
If \code{plotfn=T}, the test function and data set are plotted.
}
\section{RELEASE}{
Part of the CThresh addon to WaveThresh. Copyright Stuart Barber and Guy Nason 2004. 
}
\keyword{manip}
\author{Stuart Barber}
