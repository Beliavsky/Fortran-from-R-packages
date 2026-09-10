\name{accessD}
\alias{accessD}
\title{Get "detail" (mother wavelet) coefficients data from wavelet object}
\description{
This generic function extracts detail from various types of wavelet objects. It extracts and returns a whole resolution level of coefficients. To obtain individual packets from relevant transforms use the \link{getpacket}() series of functions. 
This function is generic. 

Particular methods exist. For objects of class:

\describe{
\item{wd}{use the \code{\link{accessD.wd}} method}
\item{wd3D}{use the \code{\link{accessD.wd3D}} method}
\item{wp}{use the \code{\link{accessD.wp}} method}
\item{wpst}{use the \code{\link{accessD.wpst}} method}
\item{wst}{use the \code{\link{accessD.wst}} method}
}

See individual method help pages for operation and examples. 
}
\usage{
accessD(\dots)
}
\arguments{
\item{\dots}{See individual help for details.}
}
\value{
A vector coefficients representing the detail coefficients for the requested resolution level. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994}
\seealso{
\code{\link{accessD.wd}}, \code{\link{accessD.wp}},\code{\link{accessD.wst}},\code{\link{accessC}}
}
\keyword{manip}
\author{G P Nason}
