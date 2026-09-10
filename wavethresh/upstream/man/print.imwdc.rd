\name{print.imwdc}
\alias{print.imwdc}
\title{Print out information about an imwdc object in readable form. }
\description{
This function prints out information about an \code{\link{imwdc.object}} in a nice human-readable form. 

Note that this function is automatically called by SPlus whenever the name of an \code{\link{imwdc.object}} is typed or whenever such an object is returned to the top level of the S interpreter. 
}
\usage{
\method{print}{imwdc}(x, ...)
}
\arguments{
\item{x}{An object of class imwdc that you wish to print out.}
\item{\dots}{This argument actually does nothing in this function! }
}
\details{
Prints out information about \code{imwdc} objects in nice readable format. 
}
\value{
The last thing this function does is call \code{\link{summary.imwdc}} so the return value is whatever is returned by this function. 
}
\section{RELEASE}{Version 2.2 Copyright Guy Nason 1994 }
\seealso{
\code{\link{imwdc.object}}, \code{\link{summary.imwdc}}. 
}
\examples{
#
# Generate an imwd object.
#
tmp <- imwd(matrix(0, nrow=32, ncol=32))
#
# Now get R to use print.imwd
#
tmp
# Class 'imwd' : Discrete Image Wavelet Transform Object:
#        ~~~~  : List with 27 components with names
#              nlevelsWT fl.dbase filter type bc date w4L4 w4L1 w4L2 w4L3
# w3L4 w3L1 w3L2 w3L3 w2L4 w2L1 w2L2 w2L3 w1L4 w1L1 w1L2 w1L3 w0L4 w0L1
# w0L2 w0L3 w0Lconstant 
#
# $ wNLx are LONG coefficient vectors !
#
# summary(.):
# ----------
# UNcompressed image wavelet decomposition structure
# Levels:  5 
# Original image was 32 x 32  pixels.
# Filter was:  Daub cmpct on least asymm N=10 
# Boundary handling:  periodic 
}
\keyword{utilities}
\author{G P Nason}
