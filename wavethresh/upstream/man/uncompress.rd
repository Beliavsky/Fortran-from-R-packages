\name{uncompress}
\alias{uncompress}
\title{Uncompress objects}
\description{
Uncompress objects.
 
This function is generic. 

Particular methods exist. For the \code{\link{imwdc.object}} class object this
generic function uses \code{\link{uncompress.imwdc}}.
There is a default uncompression method:
\code{\link{uncompress.default}} that works on vectors. 
}
\usage{
uncompress(\dots)
}
\arguments{
\item{\dots}{See individual help pages for details. }
}
\details{
See individual method help pages for operation and examples
}
\value{
A uncompressed version of the input. 
}
\section{RELEASE}{Version 2.0 Copyright Guy Nason 1993}
\seealso{
\code{\link{uncompress.default}}, \code{\link{uncompress.imwdc}}, \code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{imwdc.object}}, \code{\link{threshold.imwd}} 
}
\keyword{manip}
\keyword{utilities}
\author{G P Nason}
