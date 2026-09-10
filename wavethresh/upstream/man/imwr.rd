\name{imwr}
\alias{imwr}
\title{Inverse two-dimensional wavelet transform. }
\description{
Perform inverse two-dimensional wavelet transform using Mallat's, 1989 algorithm. 

This function is generic. 

Particular methods exist. For the \code{\link{imwd}} class object this generic function uses \code{\link{imwr.imwd}}. For the \code{imwdc} class object this generic function uses \code{\link{imwr.imwdc}}. 
}
\usage{
imwr(...)
}
\arguments{
\item{\dots}{See individual help pages for details.} 
}
\details{
See individual method help pages for operation and examples. 
}
\value{
A square matrix whose side length is a power of two that represents the inverse 2D wavelet transform of the input object x. 
}
\section{RELEASE}{Version 2 Copyright Guy Nason 1993 }
\seealso{
\code{\link{imwd}}, \code{\link{imwr.imwd}}, \code{\link{imwr.imwdc}}. 
}
\keyword{nonlinear}
\keyword{smooth}
\author{G P Nason}
