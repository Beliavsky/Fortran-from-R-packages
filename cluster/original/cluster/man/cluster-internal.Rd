\name{cluster-internal}
\alias{meanabsdev}
\title{Internal cluster functions}
\description{
  Internal cluster functions.
}
\usage{
meanabsdev(y)
}
\arguments{
  \item{y}{numeric vector.}
}
\details{
  These are not to be called by the user.

  A more \R-like version of \code{meanabsdev()} would default to
  \code{na.rm=FALSE}, as e.g., \code{\link{mean}()} or \code{\link{median}()}:
\preformatted{
  meanabsdev <- function(y, na.rm=FALSE)
       mean(abs(y - mean(y, na.rm=na.rm)), na.rm=na.rm)
}
}
\keyword{internal}
