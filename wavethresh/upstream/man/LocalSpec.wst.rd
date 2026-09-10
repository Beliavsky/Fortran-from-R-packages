\name{LocalSpec.wst}
\alias{LocalSpec.wst}
\title{Obsolete function (use ewspec)}
\usage{
\method{LocalSpec}{wst}(wst, \dots)
}
\arguments{
\item{wst}{The wst object to perform local spectral analysis on}
\item{\dots}{Other arguments to \code{\link{LocalSpec.wd}}.  }
}
\description{
This function computes a local spectra as described in Nason and Silverman
(1995). However, the function is obsolete and superceded by
\code{\link{ewspec}}.
}
\details{
Description says it all.

However, this function converts the \code{\link{wst.object}} object to
a nondecimated \code{\link{wd.object}} and then calls
\code{\link{LocalSpec.wd}}.
}
\value{
Same value as \code{\link{LocalSpec.wd}}.
}
\seealso{\code{\link{ewspec}}}
\author{G P Nason}
\keyword{ts}
