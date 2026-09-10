\name{mpostfilter}
\alias{mpostfilter}
\title{Multiwavelet postfilter}
\usage{
mpostfilter(C, prefilter.type, filter.type, 
    nphi, npsi, ndecim, nlevels, verbose = FALSE)

}
\arguments{
\item{C}{The multivariate sequence you wish to turn back into a univariate
	one using the inverse of an earlier prefilter operation.}
\item{prefilter.type}{Controls the type of prefilter (see Tim Downie's 
	PhD thesis, or references therein. Types include \code{Minimal},
	\code{Identity}, \code{Repeat}, \code{Interp}, \code{default},
	\code{Xia}, \code{Roach1}, \code{Roach3}, \code{Donovan3} or
	\code{Linear}}
\item{filter.type}{The type of multiwavelet: can be \code{Geronimo}
	or \code{Donovan3}}
\item{nphi}{The number of father wavelets in the system}
\item{npsi}{The number of mother wavelets in the system}
\item{ndecim}{The ndecim parameter (not apparently used here)}
\item{nlevels}{The number of levels in the multiwavelet transform}
\item{verbose}{If TRUE then informative messages are printed as the
	function progresses}
}
\description{
A multiwavelet postfilter turns a multivariate sequence into a univariate
sequence. As such, the postfilter is used
on the inverse transform, it is the inverse of an earlier used prefilter.

Not intended for direct user use.
}
\details{
Description says all
}
\value{
The appropriate postfiltered data.
}
\seealso{\code{\link{mprefilter}},\code{\link{mwd}}}
\author{Tim Downie}
\keyword{math}
