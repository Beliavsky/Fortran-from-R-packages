\name{wd3D.object}
\alias{wd3D.object}
\title{Three-dimensional wavelet object}
\description{
These are objects of classes 

wd3D

They contain the 3D discrete wavelet transform of a 3D array (with each dimension being the same dyadic size). 
}
\details{
To retain your sanity the wavelet coefficients at any resolution level in directions, GGG, GGH, GHG, GHH, HGG, HGH, HHG should be extracted by the \code{\link{accessD}}() function and inserted using the \code{\link{putD}} function rather than by the \code{$} operator.
}
\value{
The following components must be included in a legitimate `wd' object. 

\item{a}{a three-dimensional array containing the 3D discrete wavelet coefficients. The coefficients are stored in a pyramid structure for efficiency.}
\item{nlevelsWT}{The number of levels in the pyramidal decomposition that produces the coefficients. If you raise 2 to the power of nlevels you get the number of data points used in each dimension of the decomposition.}
\item{filter.number}{the number of the wavelet family that did the DWT.} 
\item{family}{the family of wavelets that did the DWT.} 
\item{date}{the date that the transform was computed.}
}
\section{generation}{
This class of objects is returned from the wd3D function to represent a three-dimensional DWT of a 3D array. Other functions return an object of class wd3D.
}
\section{methods}{
The wd3D class of objects has methods for the following generic functions: \code{\link{accessD}}, \code{\link{print}}, \code{\link{putD}}, \code{\link{summary}}, \code{\link{threshold}}.
}
\section{release}{Version 3.9.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{wd3D}}, \code{\link{accessD.wd3D}}, \code{\link{print.wd3D}}, \code{\link{putD.wd3D}}, \code{\link{putDwd3Dcheck}}, \code{\link{summary.wd3D}}, \code{\link{threshold.wd3D}}, \code{\link{wr3D}}. 
}
\keyword{classes}
\author{G P Nason}
