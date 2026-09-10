\name{putC}
\alias{putC}
\title{Put smoothed data (father wavelet) coefficients into wavelet structure}
\description{
This generic function inserts smooths into various types of wavelet objects. 

This function is generic. 

Particular methods exist. For objects of class: 

\describe{
\item{wd}{use the \code{\link{putC.wd}} method.}
\item{wp}{use the \code{\link{putC.wp}} method.} 
\item{wst}{use the \code{\link{putC.wst}} method.}
}
 
See individual method help pages for operation and examples. 

See \code{\link{accessC}} if you wish to \emph{extract} father wavelet coefficients. See \code{\link{putD}} if you wish to insert \emph{mother} wavelet coefficients
}
\usage{
putC(...)
}
\arguments{
\item{\dots}{See individual help pages for details.} 
}
\value{
A wavelet object of the same class as x with the new father wavelet coefficients inserted. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{putC.wd}}, \code{\link{putC.wp}}, \code{\link{putC.wst}}, \code{\link{accessC}}, \code{\link{putD}}. 
}
\keyword{manip}
\author{G P Nason}
