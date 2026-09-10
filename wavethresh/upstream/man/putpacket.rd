\name{putpacket}
\alias{putpacket}
\title{Insert a packet of coefficients into a wavelet object. }
\description{
This generic function inserts packets of coefficients into various types of wavelet objects. 

This function is generic. 

Particular methods exist. For objects of class: 

\describe{
\item{wp}{use the \code{\link{putpacket.wp}} method.}
\item{wst}{use the \code{\link{putpacket.wst}} method.} 
\item{wst2D}{use the \code{\link{putpacket.wst2D}} method.}
}
 
See individual method help pages for operation and examples. 

Use the \code{\link{putC}} and \code{\link{putD}} function to insert whole resolution levels of coefficients simultaneously. 
}
\usage{
putpacket(\dots)
}
\arguments{
\item{\dots}{See individual help pages for details.}
}
\value{
A wavelet object of the same class as \code{x} the input object. The returned wavelet object is the same as the input except that the appropriate packet of coefficients supplied is replaced. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{putpacket.wp}}, \code{\link{putpacket.wst}}, \code{\link{putpacket.wst2D}}, \code{\link{putD}}, \code{\link{putC}}, \code{\link{wp.object}}, \code{\link{wst.object}}, \code{\link{wst2D.object}}. 
}
\keyword{manip}
\author{G P Nason}
