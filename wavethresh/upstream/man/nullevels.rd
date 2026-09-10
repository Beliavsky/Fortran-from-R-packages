\name{nullevels}
\alias{nullevels}
\title{Set whole resolution levels of coefficients equal to zero.}
\description{
Generic function which sets whole resolution levels of coefficients equal to zero. 

Particular methods exist. For objects of class: 

\describe{
\item{imwd}{use the \code{\link{nullevels.imwd}} method.} 
\item{wd}{use the \code{\link{nullevels.wd}} method. }
\item{wst}{use the \code{\link{nullevels.wst}} method.} 
}

See individual method help pages for operation and examples. 
}
\usage{
nullevels(\dots)
}
\arguments{
\item{\dots}{See individual help pages for details.}
}
\value{
An object of the same class as x but with the specified levels set to zero. 
}
\section{RELEASE}{Version 3.8.1 Copyright Guy Nason 1997 }
\seealso{
\code{\link{nullevels.imwd}} \code{\link{nullevels.wd}} \code{\link{nullevels.wst}} \code{\link{wd.object}}, \code{\link{wd}} \code{\link{wst.object}} \code{\link{wst}} 
}
\keyword{manip}
\author{G P Nason}
