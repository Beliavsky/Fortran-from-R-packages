\name{print.wd3D}
\alias{print.wd3D}
\title{Print out information about an wd3D object in a readable form.}           
\description{
This function prints out information about an \code{\link{wd3D.object}} in a readable form. 

Note that this function is automatically called by SPlus whenever the name of an \code{\link{wd3D.object}} is typed or whenever such an object is returned to the top level of the S interpreter

}
\usage{
\method{print}{wd3D}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link{wd3D}} that you wish to print out.}
\item{\dots}{This argument actually does nothing in this function!}
}
\details{
Prints out information about \code{\link{wd3D}} objects in nice readable format. 
}
\value{
The last thing this function does is call \code{\link{summary.wd3D}} so the return value is whatever is returned by this function.} 

\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{accessD.wd3D}}, \code{\link{print.wd3D}}, \code{\link{putD.wd3D}}, \code{\link{putDwd3Dcheck}}, \code{\link{summary.wd3D}}, \code{\link{threshold.wd3D}}, \code{\link{wd3D}}, \code{\link{wd3D.object}}, \code{\link{wr3D}}. 

}
\examples{
#
# Generate an wd3D object.
#
tmp <- wd3D(array(rnorm(512), dim=c(8,8,8)))
#
# Now get R to use print.wd
#
tmp
#Class 'wd3d' : 3D DWT Object:
#       ~~~~  : List with 5 components with names
#              a filter.number family date nlevelsWT 
#
#$ a is the wavelet coefficient array
#Dimension of a is [1] 8 8 8
#
#Created on : Wed Oct 20 17:24:15 BST 1999 
#
#summary(.):
#----------
#Levels:  3 
#Filter number was:  10 
#Filter family was:  DaubLeAsymm 
#Date:  Wed Oct 20 17:24:15 BST 1999 
}
\keyword{manip}
\author{G P Nason}
