\name{levarr}
\alias{levarr}
\title{Subsidiary routine that generates a particular permutation}
\usage{
levarr(v, levstodo)
}
\arguments{
\item{v}{the vector to permute}
\item{levstodo}{the number of levels associated with the current
	level in the object you wish to permute}
}
\description{
Not intended for casual user use. This function is used to provide
the partition to reorder \code{\link{wst.object}} into \code{\link{wd.object}}
(nondecimated time ordered) objects.
}
\details{
Description says all
}
\value{
A permutation of the \code{v} vector according to the number of levels
that need handling
}
\seealso{\code{\link{getarrvec}}, \code{\link{convert.wd}}, \code{\link{convert.wst}}}
\examples{
levarr(1:4, 3)
# [1] 1 3 2 4
}
\author{G P Nason}
\keyword{manip}
