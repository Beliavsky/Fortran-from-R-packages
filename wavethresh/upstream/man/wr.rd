\name{wr}
\alias{wr}
\title{Wavelet reconstruction (inverse DWT).}
\description{
Performs inverse discrete wavelet transform. 

This function is generic. 

Particular methods exist. For the \code{\link{wd}} class object this generic function uses \code{\link{wr.wd}}. 

}
\usage{
wr(...)
}
\arguments{
\item{\dots}{See individual help pages for details.}
}
\details{
See individual method help pages for operation and examples. 
}
\value{
Usually the wavelet reconstruction of x. Although the return value varies with the precise method used. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wd}}, \code{\link{wd.object}}, \code{\link{wr.wd}} 
}
\keyword{manip}
\author{G P Nason}
