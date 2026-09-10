\name{nlevelsWT}
\alias{nlevelsWT}
\title{Returns number of scale (resolution) levels.}
\description{
Returns the number of scales (or resolutions) in various wavelet objects and for some objects returns the number of scales that would result if processed by a wavelet routine. 

This function is generic. 

One methods exists at present as most wavelet objects store the number of levels as the \code{nlevelsWT} component. The method that exists is\code{\link{nlevelsWT.default}} 
}
\usage{
nlevelsWT(...)
}
\arguments{
\item{\dots}{See individual help pages for details.} 
}
\details{
See individual method help pages for operation and examples. 
}
\value{
An integer representing the number of levels associated with the object. 
}
\section{RELEASE}{Version 3.6.0 Copyright Guy Nason 1995 }
\seealso{
\code{\link{nlevelsWT.default}}
}
\keyword{arith}
\author{G P Nason}
