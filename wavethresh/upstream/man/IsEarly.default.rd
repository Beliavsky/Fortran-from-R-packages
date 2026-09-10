\name{IsEarly.default}
\alias{IsEarly.default}
\title{Detects whether object is from an earlier version of WaveThresh}
\usage{
\method{IsEarly}{default}(x)
}
\arguments{
\item{x}{Object to discern}
}
\description{
Detects whether object is from an earlier version of WaveThresh.
}
\details{
The default method always returns FALSE, i.e. unless the object is of
a specific type handled by a particular method then it won't be from an
earlier version.
}
\value{
Always FALSE for the generic
}
\seealso{\code{\link{IsEarly}}}
\author{G P Nason}
\keyword{error}
