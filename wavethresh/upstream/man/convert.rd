\name{convert}
\alias{convert}
\title{Convert one type of wavelet object into another. }
\description{
Convert one type of wavelet object into another. 

This function is generic. 

Particular methods exist: 

\code{\link{convert.wd}} is used to convert non-decimated \code{\link{wd}} objects into \code{\link{wst}} objects. 
\code{\link{convert.wst}} is used to convert \code{\link{wst}} objects into non-decimated \code{\link{wd}} objects. 
}
\usage{
convert(...)
}
\arguments{
\item{\dots}{See individual help pages for details.}
}
\details{
See individual method help pages for operation and examples. 
}
\value{
An object containing the converted representation. 
}
\section{RELEASE}{Version 3.6.0 Copyright Guy Nason 1995 }
\seealso{
\code{\link{convert.wd}}, \code{\link{convert.wst}}, \code{\link{wst}}, \code{\link{wst.object}}, \code{\link{wst}}, \code{\link{wst.object}}. 
}

\keyword{manip}
\author{G P Nason}
