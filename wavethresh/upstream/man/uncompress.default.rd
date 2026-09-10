\name{uncompress.default}
\alias{uncompress.default}
\title{Undo zero run-length encoding for a vector.}
\usage{
\method{uncompress}{default}(v, verbose=FALSE, ...)
}
\arguments{
\item{v}{The object to uncompress}
\item{verbose}{Print an informative message whilst executing}
\item{\dots}{Other arguments}
}
\description{
This function inverts the action carried out by the
\code{\link{compress.default}} function.
}
\details{
The inverse of \code{\link{compress.default}}
}
\value{
The uncompressed, reinstated, vector.
}
\seealso{\code{\link{compress.default}}, \code{\link{uncompress}}}
\examples{
uncompress(compress(c(1, rep(0,99), 1)))
#[1] 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
}
\author{G P Nason}
\keyword{manip}
