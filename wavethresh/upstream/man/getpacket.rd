\name{getpacket}
\alias{getpacket}
\title{Get a packet of coefficients from a wavelet object}
\description{
This generic function extracts packets of coefficients from various types of wavelet objects. 

This function is generic. 

Particular methods exist. For objects of class: 

\describe{
\item{wp}{use the \code{\link{getpacket.wp}} method.} 
\item{wst}{use the \code{\link{getpacket.wst}} method.} 
\item{wpst}{use the \code{\link{getpacket.wpst}} method.}
}
 
See individual method help pages for operation and examples. 

Use the \code{\link{accessC}} and \code{\link{accessD}} function to extract whole resolution levels of coefficients simultaneously. 
}
\usage{
getpacket(...)
}
\arguments{
\item{\dots}{See individual help pages for details.} 
}
\value{
The packet of coefficients requested. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{getpacket.wp}}, \code{\link{getpacket.wst}}, \code{\link{getpacket.wpst}}, \code{\link{accessD}}, \code{\link{accessC}}. 
}
\keyword{manip}
\author{G P Nason}
