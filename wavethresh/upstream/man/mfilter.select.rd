\name{mfilter.select}
\alias{mfilter.select}
\title{Provide filter coefficients for multiple wavelets.}
\description{
This function returns the filter coefficients necessary for doing a discrete multiple wavelet transform (and its inverse). 
}
\usage{
mfilter.select(type = "Geronimo")
}
\arguments{
\item{type}{The name for the multiple wavelet basis. The two possible types are "Geronimo" and "Donovan3"}. 
}
\details{
This function supplies the multiple wavelet filter coefficients required by the \code{\link{mwd}} function.
 
A multiple wavelet filter is somewhat different from a single wavelet filter. Firstly the filters are made up of matrices not single coefficients. Secondly there is no simple expression for the high pass coefficients G in terms of the low pass coefficients H, so both sets of coefficients must be specified. Note also that the transpose of the filter coefficients are used in the inverse transform, an unnecessary detail with scalar coefficients. There are two filters available at the moment. Geronimo is the default, and is recommended as it has been checked thoroughly. Donovan3 uses three orthogonal wavelets described in Donovan et al. but this coding has had little testing. 

See Donovan, Geronimo and Hardin, 1996 and Geronimo, Hardin and Massopust, 1994. 

This function fulfils the same purpose as the \code{\link{filter.select}} function does for the standard DWT \code{\link{wd}}. 

}
\value{
A list is returned with the following eight components which describe the filter:
 
\item{type}{The multiple wavelet basis type string.} 
\item{H}{A vector containing the low pass filter coefficients.}
\item{G}{A vector containing the high pass pass filter coefficients.} 
\item{name}{A character string containing the full name of the filter.}
\item{nphi}{The number of scaling functions in the multiple wavelet basis.} 
\item{npsi}{The number of wavelet functions in the multiple wavelet basis.}
\item{NH}{The number of matrix coefficients in the filter. This is different from length(H).}
\item{ndecim}{The decimation factor. I.e. the scale ratio between two successive resolution levels.}
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6)}
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mwd.object}}, \code{\link{mwd}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#This function is currently used by `mwr' and `mwd' in decomposing and
#reconstructing, however you can view the coefficients.
#
# look at the filter coefficients for Geronimo multiwavelet
#
mfilter.select()
#$type:
#[1] "Geronimo"
#
#$name:
#[1] "Geronimo Multiwavelets"
#
#$nphi:
#[1] 2
#
#$npsi:
#[1] 2
#
#$NH:
#[1] 4
#
#$ndecim:
#[1] 2
#$H:
# [1]  0.4242641  0.8000000 -0.0500000 -0.2121320  0.4242641  0.0000000
# [7]  0.4500000  0.7071068  0.0000000  0.0000000  0.4500000 -0.2121320
#[13]  0.0000000  0.0000000 -0.0500000  0.0000000
#
#$G:
# [1] -0.05000000 -0.21213203  0.07071068  0.30000000  0.45000000 -0.70710678
#
# [7] -0.63639610  0.00000000  0.45000000 -0.21213203  0.63639610 -0.30000000
#[13] -0.05000000  0.00000000 -0.07071068  0.00000000

}
\keyword{manip}
\author{Tim Downie}
