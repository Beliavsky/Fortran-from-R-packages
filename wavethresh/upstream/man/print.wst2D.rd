\name{print.wst2D}
\alias{print.wst2D}
\title{Print out information about an wst2d object in a readable form.}           
\description{
This function prints out information about an \code{\link{wst2D.object}} in a nice human- readable form. 

Note that this function is automatically called by SPlus whenever the name of an \code{\link{wst2D.object}} is typed or whenever such an object is returned to the top level of the S interpreter

}
\usage{
\method{print}{wst2D}(x, ...)
}
\arguments{
\item{x}{An object of class \code{\link{wst2D}} that you wish to print out.}
\item{\dots}{This argument actually does nothing in this function!}
}
\details{
Prints out information about \code{\link{wst2D}} objects in nice readable format. 
}
\value{
The last thing this function does is call \code{\link{summary.wst2D}} so the return value is whatever is returned by this function.} 

\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1998 }
\seealso{
\code{\link{wst2D.object}}, \code{\link{summary.wst2D}}. 
}
\examples{
#
# This examples uses the uawst2D object created in the EXAMPLES
# section of the help page for wst2D 
#
#uawst2D
#Class 'wst2D' : 2D Stationary Wavelet Transform Object:
#       ~~~~~  : List with 4 components with names
#              wst2D nlevelsWT filter date 
#
#$wst2D is the coefficient array
#
#Created on : Fri Nov  5 18:06:17 GMT 1999 
#
#summary(.):
#----------
#Levels:  8 
#Length of original:  256 x 256 
#Filter was:  Daub cmpct on least asymm N=10 
#Date:  Fri Nov  5 18:06:17 GMT 1999 
}
\keyword{utilities}
\author{G P Nason}
