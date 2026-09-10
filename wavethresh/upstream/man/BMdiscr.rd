\name{BMdiscr}
\alias{BMdiscr}
\title{Subsidiary routine for makewpstDO function}
\description{
Function actually performs discrimination on reduced variable set supplied to
it from \code{\link{Best1DCols}} function.
}
\usage{
BMdiscr(BP)
}
\arguments{
\item{BP}{An list of the same format as returned by \code{\link{Best1DCols}}}
}
\details{
Not intended for direct user use
}
\value{
Returns a list of objects: essentially the input argument \code{BP}
and the return value from a call to the \code{lda} function which
performs the discrimination operation.
}
\seealso{\code{\link{Best1DCols}},\code{\link{makewpstDO}}}
\author{G P Nason}
\keyword{ts}
\keyword{multivariate}
