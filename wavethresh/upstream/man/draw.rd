\name{draw}
\alias{draw}
\title{Draw wavelets or scaling functions.}
\description{
Draws the mother wavelet or scaling function associated with an object. 

This function is generic. 

Particular methods exist. The following functions are used for the following objects: 

\describe{
\item{imwd.object}{the \code{\link{draw.imwd}} function is used.}
\item{imwdc.object}{the \code{\link{draw.imwdc}} function is used.}
\item{wd.object}{the \code{\link{draw.wd}} function is used.}
\item{wp.object}{the \code{\link{draw.wp}} function is used.}
\item{wst.object}{the \code{\link{draw.wst}} function is used.}
}

All of the above method functions use the \code{\link{draw.default}} function which is the function which actually does the drawing. 
}
\usage{
draw(...)
}
\arguments{
\item{\dots}{methods may have additional arguments}
}
\details{
See individual method help pages for operation and examples.
}
\value{
If the \code{plot.it} argument is supplied then the draw functions tend to return the coordinates of what they were meant to draw and don't actually draw anything. 
}
\section{RELEASE}{Version 2 Copyright Guy Nason 1993}
\seealso{
\code{\link{draw.default}}, \code{\link{draw.imwd}}, \code{\link{draw.imwdc}}, \code{\link{draw.wd}}, \code{\link{draw.wp}}, \code{\link{draw.wst}}, \code{\link{imwd.object}}, \code{\link{imwdc.object}}, \code{\link{wd.object}}, \code{\link{wp.object}}, \code{\link{wst.object}}. 

}
\keyword{hplot}
\author{G P Nason}
