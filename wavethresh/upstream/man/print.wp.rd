\name{print.wp}
\alias{print.wp}
\title{Print out information about an wd object in readable form. }
\description{
This function prints out information about an \code{\link{wp.object}} in a nice human-readable form. 

Note that this function is automatically called by SPlus whenever the name of an \code{\link{wp.object}} is typed or whenever such an object is returned to the top level of the S interpreter

}
\usage{
\method{print}{wp}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link{wp}} that you wish to print out.}
\item{\dots}{This argument actually does nothing in this function!}
}
\details{
Prints out information about \code{\link{wp}} objects in nice readable format. 
}
\value{
The last thing this function does is call \code{\link{summary.wp}} so the return value is whatever is returned by this function.} 

\section{RELEASE}{Version 3.0 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wp.object}}, \code{\link{summary.wp}}. 
}
\examples{
#
# Generate an wp object.
#
tmp <- wp(rnorm(32))
#
# Now get Splus to use print.wp
#
tmp
#
# Now get Splus to use print.wp
#
# tmp
# Class 'wp' : Wavelet Packet Object:
#        ~~  : List with 4 components with names
#               wp nlevelsWT filter date 
# 
# $wp is the wavelet packet matrix
# 
# Created on : Fri Oct 23 19:59:01 1998 
# 
# summary(.):
# ----------
# Levels:  5 
# Length of original:  32 
# Filter was:  Daub cmpct on least asymm N=10 
}
\keyword{manip}
\author{G P Nason}
