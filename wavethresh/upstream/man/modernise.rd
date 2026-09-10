\name{modernise}
\alias{modernise}
\title{Generic function to upgrade a V2 WaveThresh object to V4}
\usage{
modernise(...)
}
\arguments{
\item{\dots}{Other objects}
}
\description{
Not really used in practice. The function \code{\link{IsEarly}}
can be used to tell if an object comes from an earlier version of
wavethresh. Note that the earlier version only has a \code{\link{wd.object}}
class object so there is only a method for that.
}
\details{
Description says all
}
\value{
A modernised version of the object.
}
\seealso{\code{\link{IsEarly}},\code{\link{modernise.wd}}}
\author{G P Nason}
\keyword{manip}
