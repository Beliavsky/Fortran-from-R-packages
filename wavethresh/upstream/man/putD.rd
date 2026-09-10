\name{putD}
\alias{putD}
\title{Put mother wavelet coefficients into wavelet structure}
\description{
This generic function inserts smooths into various types of wavelet objects. 

This function is generic. 

Particular methods exist. For objects of class: 

\describe{
\item{wd}{use the \code{\link{putD.wd}} method.}
\item{wp}{use the \code{\link{putD.wp}} method. }
\item{wst}{use the \code{\link{putD.wst}} method.}
}

 
See individual method help pages for operation and examples. 

See \code{\link{accessD}} if you wish to \emph{extract} mother wavelet coefficients. See \code{\link{putC}} if you wish to insert \emph{father} wavelet coefficients. 
}
\usage{
putD(\dots)
}
\arguments{
\item{\dots}{See individual help pages for details.} 
}
\value{
A wavelet object of the same class as \code{x} with the new mother wavelet coefficients inserted. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{putD.wd}}, \code{\link{putD.wp}}, \code{\link{putD.wst}}, \code{\link{accessD}}, \code{\link{putC}}. 
}
\keyword{manip}
\author{G P Nason}
