\name{mprefilter}
\alias{mprefilter}
\title{Multiwavelet prefilter}
\usage{
mprefilter(data, prefilter.type, filter.type, nlevels, nvecs.c,
    nphi, npsi, ndecim, verbose = FALSE)

}
\arguments{
\item{data}{The univariate sequence that you wish to turn into a multivariate one}
\item{prefilter.type}{Controls the type of prefilter (see Tim Downie's 
	PhD thesis, or references therein. Types include \code{Minimal},
	\code{Identity}, \code{Repeat}, \code{Interp}, \code{default},
	\code{Xia}, \code{Roach1}, \code{Roach3}, \code{Donovan3} or
	\code{Linear}}
\item{filter.type}{The type of multiwavelet: can be \code{Geronimo}
	or \code{Donovan3}}
\item{nlevels}{The number of levels in the multiwavelet transform}
\item{nvecs.c}{Parameter obtained from the mfirst.last function
	related to the particular filters}
\item{nphi}{The number of father wavelets in the system}
\item{npsi}{The number of mother wavelets in the system}
\item{ndecim}{The ndecim parameter (not apparently used here)}
\item{verbose}{If TRUE then informative messages are printed as the
	function progresses}
}
\description{
A multiwavelet prefilter turns a univariate sequence into a bivariate
(in this case) sequence suitable for processing by a multiwavelet
transform, such as \code{\link{mwd}}. As such, the prefilter is used
on the forward transform.

Not intended for direct user use.
}
\details{
Description says all
}
\value{
The appropriate prefiltered data.
}
\seealso{\code{\link{mpostfilter}},\code{\link{mwd}}}
\author{Tim Downie}
\keyword{math}
