\name{AvBasis}
\alias{AvBasis}
\title{Basis averaging ("inversion")}
\description{
Average of whole collection of basis functions. 

This function is generic. 

Particular methods exist. For the \code{\link{wst}} class object this generic function uses  \code{\link{AvBasis.wst}}. In the future we hope to add methods for  \code{\link{wp}} and  \code{\link{wpst}} class objects.
}
\usage{
AvBasis(...)
}
\arguments{
\item{\dots}{See individual help pages for details}
}
\details{
See individual method help pages for operation and examples. 
}
\value{
A vector containing the average of the representation over all bases. 
}
\section{RELEASE}{
Version 3.6.0 Copyright Guy Nason 1995 
}
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{AvBasis.wst}}
}
\keyword{manip} 
\author{G P Nason}
