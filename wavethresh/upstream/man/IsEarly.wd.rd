\name{IsEarly.wd}
\alias{IsEarly.wd}
\title{Function to detect whether a wd object is from WaveThresh2 or not}
\usage{
\method{IsEarly}{wd}(x)
}
\arguments{
\item{x}{The wd object that you are trying to check}
}
\description{
Function to detect whether a wd object is from WaveThresh2 or not.
}
\details{
The function merely looks to see whether the wd object has a component
called date. If it does not then it is from version 2. This routine
is legacy and not very important anymore.
}
\value{
Returns TRUE if from an earlier version of WaveThresh (v2), returns
FALSE if not. 
}
\seealso{\code{\link{IsEarly}}}
\author{G P Nason}
\keyword{misc}
