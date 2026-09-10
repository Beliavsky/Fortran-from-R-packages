\name{summary.mwd}
\alias{summary.mwd}
\title{Use summary() on a mwd object. }
\description{
This function prints out more information about an \code{\link{mwd.object}} in a nice human-readable form. 
}
\usage{
\method{summary}{mwd}(object, ...)
}
\arguments{
\item{object}{An object of class \code{\link{mwd}} that you wish to print out more information.}
\item{...}{Any other arguments.}

}
\value{
Nothing of any particular interest.
}
\note{Prints out information about \code{\link{mwd}} objects in nice readable format.} 
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6) }
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Generate an mwd object.
#
tmp <- mwd(rnorm(32))
#
# Now get Splus to use summary.mwd
#
summary(tmp)
# Length of original:  32 
# Levels:  4 
# Filter was:  Geronimo Multiwavelets 
# Scaling fns:  2 
# Wavelet fns:  2 
# Prefilter:  default 
# Scaling factor:  2 
# Boundary handling:  periodic 
# Transform type:  wavelet 
# Date:  Tue Nov 16 13:55:26 GMT 1999 
}
\keyword{nonlinear}
\keyword{smooth}
\author{Tim Downie}
