\name{accessc}
\alias{accessc}
\title{Get variance information from irregularly spaced wavelet decomposition object. }
\description{
This function gets information from the c component of an \code{\link{irregwd.objects}} an irregularly spaced wavelet decomposition object. 

Note that this function is \emph{not} the same as \code{\link{accessC}} which obtains father wavelet coefficients from an \code{\link{wd}} class object. 
}
\usage{
accessc(irregwd.structure, level, boundary=FALSE)
}
\arguments{
\item{irregwd.structure}{Irregular wavelet decomposition object from which you wish to extract parts of the \code{c} component from.}
\item{level}{The level that you wish to extract. This value ranges from 0 to the \code{\link{nlevelsWT}}(irregwd.structure)-1.}
\item{boundary}{If this argument is T then all of the boundary correction values will be returned as well (note: the length of the returned vector may not be a power of 2). If boundary is false, then just the coefficients will be returned. If the decomposition (or reconstruction) was done with periodic boundary conditions then this option has no effect.}
}
\details{
The \code{\link{irregwd}} function produces a irregular wavelet decomposition (reconstruction) structure. 

The \code{c} component is stored in a similar way to the C and D vectors which store the father and mother wavelet coefficients respectively. Hence to access the information the accessc function plays a similar role to \code{\link{accessC}} and \code{\link{accessD}} functions. 
}
\value{
A vector of the extracted data. 
}
\section{RELEASE}{Version 3.9.4 Code Copyright Arne Kovac 1997. Help Copyright Guy Nason 2004. }
\seealso{
\code{\link{irregwd}}, \code{\link{irregwd.objects}}, \code{\link{threshold.irregwd}},\code{\link{makegrid}}, \code{\link{plot.irregwd}}}
\examples{
#
# Most users will not need to use this function. However, see the main
# examples for the irregular wavelet denoising in the examples for
# makegrid.
#
}
\keyword{manip}
\author{G P Nason}
