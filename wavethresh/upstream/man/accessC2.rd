\name{accessC}
\alias{accessC}
\title{Get "detail" (mother wavelet) coefficients data from wavelet object}
\description{
This generic function extracts detail from various types of wavelet objects. It extracts and returns a whole resolution level of coefficients. To obtain individual packets from relevant transforms use the \link{getpacket}() series of functions. 
This function is generic. 

Particular methods exist. For objects of class:

\describe{
\item{wd}{use the \code{\link{accessC.wd}} method}
%\item{wd3D}{use the \code{\link{accessD.wd3D}} method}
\item{wp}{use the \code{\link{accessC.wp}} method}
%\item{wpst}{use the \code{\link{accessC.wpst}} method}
\item{wst}{use the \code{\link{accessC.wst}} method}
}

See individual method help pages for operation and examples. 
}
\usage{
accessC(\dots)
}
\arguments{
\item{\dots}{See individual help for details.}
}
\value{
A vector coefficients representing the detail coefficients for the requested resolution level. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994}
\seealso{
\code{\link{accessC.wd}}, \code{\link{accessC.wp}},\code{\link{accessC.wst}},\code{\link{accessD}}
}
\keyword{manip}
\author{G P Nason}
