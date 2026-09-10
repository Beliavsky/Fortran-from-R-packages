\name{print.wd}
\alias{print.wd}
\title{Print out information about an wd object in readable form. }
\description{
This function prints out information about an \code{\link{wd.object}} in a nice human-readable form. 

Note that this function is automatically called by SPlus whenever the name of an \code{\link{wd.object}} is typed or whenever such an object is returned to the top level of the S interpreter

}
\usage{
\method{print}{wd}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link{wd}} that you wish to print out.}
\item{\dots}{This argument actually does nothing in this function!}
}
\details{
Prints out information about \code{\link{wd}} objects in nice readable format. 
}
\value{
The last thing this function does is call \code{\link{summary.wd}} so the return value is whatever is returned by this function.} 

\section{RELEASE}{Version 3.0 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wd.object}}, \code{\link{summary.wd}}. 
}
\examples{
#
# Generate an wd object.
#
tmp <- wd(rnorm(32))
#
# Now get R to use print.wd
#
tmp
# Class 'wd' : Discrete Wavelet Transform Object:
#        ~~  : List with 8 components with names
#               C D nlevelsWT fl.dbase filter type bc date 
#
# $ C and $ D are LONG coefficient vectors !
#
# Created on : Fri Oct 23 19:56:00 1998 
# Type of decomposition:  wavelet 
# 
# summary(.):
# ----------
# Levels:  5 
# Length of original:  32 
# Filter was:  Daub cmpct on least asymm N=10 
# Boundary handling:  periodic 
# Transform type:  wavelet 
# Date:  Fri Oct 23 19:56:00 1998 
#
#
}
\keyword{print}
\author{G P Nason}
