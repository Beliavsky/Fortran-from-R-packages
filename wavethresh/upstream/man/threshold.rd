\name{threshold}
\alias{threshold}
\title{Threshold coefficients}
\description{
Modify coefficients by thresholding or shrinkage. 

This function is generic. 

Particular methods exist for the following objects: 

\describe{
\item{wd object}{the \code{\link{threshold.wd}} function is used;}
\item{imwd object}{the \code{\link{threshold.imwd}} function is used;}
\item{imwdc object}{the \code{\link{threshold.imwdc}} function is used;}
\item{irregwd object}{the \code{\link{threshold.irregwd}} function is used;}
\item{wd3D object}{the \code{\link{threshold.wd3D}} function is used;} 
\item{wp object}{the \code{\link{threshold.wp}} function is used;}
\item{wst object}{the \code{\link{threshold.wst}} function is used.} 
}
}
\usage{
threshold(...)
}
\arguments{
\item{\dots}{See individual help pages for details.}
}
\details{
See individual method help pages for operation and examples. 
}
\value{
Usually a copy of the input object but containing thresholded or shrunk coefficients. 
}
\section{RELEASE}{Version 2 Copyright Guy Nason 1993 }
\seealso{
\code{\link{imwd.object}}, \code{\link{imwdc.object}}, \code{irregwd object}, \code{\link{threshold.imwd}}, \code{\link{threshold.imwdc}}, \code{\link{threshold.irregwd}}, \code{\link{threshold.wd}}, \code{\link{threshold.wd3D}}, \code{\link{threshold.wp}}, \code{\link{threshold.wst}} \code{\link{wd.object}}, \code{\link{wd3D.object}}, \code{\link{wp.object}}, \code{\link{wst.object}}. 
}
\keyword{manip}
\author{G P Nason}
