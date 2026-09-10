\name{compress}
\alias{compress}
\title{Compress objects}
\description{
Compress objects.
 
This function is generic. 

Particular methods exist. For the \code{\link{imwd}} class object this generic function uses \code{\link{compress.imwd}}. There is a default compression method: \code{\link{compress.default}} that works on vectors. 
}
\usage{
compress(\dots)
}
\arguments{
\item{\dots}{See individual help pages for details. }
}
\details{
See individual method help pages for operation and examples
}
\value{
A compressed version of the input. 
}
\section{RELEASE}{Version 2.0 Copyright Guy Nason 1993}
\seealso{
\code{\link{compress.default}}, \code{\link{compress.imwd}}, \code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{imwdc.object}}, \code{\link{threshold.imwd}} 
}
\keyword{manip}
\keyword{utilities}

\author{G P Nason}
