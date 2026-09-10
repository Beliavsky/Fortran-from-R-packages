\name{print.mwd}
\alias{print.mwd}
\title{Use print() on a mwd object.}
\description{
This function prints out information about an \code{\link{mwd.object}} in a nice human-readable form. 

Note that this function is automatically called by SPlus whenever the name of an \code{\link{mwd.object}} is typed or whenever such an object is returned to the top level of the S interpreter. 

}
\usage{
\method{print}{mwd}(x, ...)
}
\arguments{
\item{x}{An object of class mwd that you wish to print out.} 
\item{\dots}{This argument actually does nothing in this function!}
}
\details{
Prints out information about \code{\link{mwd}} objects in nice readable format. 
}
\value{
The last thing this function does is call \code{\link{summary.mwd}} so the return value is whatever is returned by this function. 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6) }
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}},\code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Generate an mwd object.
#
tmp <- mwd(rnorm(32))
#
# Now get Splus to use print.mwd
#
tmp
# Class 'mwd' : Discrete Multiple Wavelet Transform Object:
#	 ~~~  : List with 10 components with names
#		 C D nlevelsWT ndata filter fl.dbase type bc prefilter date 
# 
# $ C and $ D are LONG coefficient vectors !
# 
# Created on : Tue Nov 16 13:16:07 GMT 1999 
# Type of decomposition:  wavelet 
# 
# summary:
# ----------
# Length of original:  32 
# Levels:  4 
# Filter was:  Geronimo Multiwavelets 
# Scaling fns:  2 
# Wavelet fns:  2 
# Prefilter:  default 
# Scaling factor:  2 
# Boundary handling:  periodic 
# Transform type:  wavelet 
# Date:  Tue Nov 16 13:16:07 GMT 1999 
}
\keyword{utilities}
\author{G P Nason}
