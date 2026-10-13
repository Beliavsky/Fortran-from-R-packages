\name{major}
\alias{major}
\title{Majorization}
\description{
 tests whether a vector \code{x} majorizes another vector \code{y}}
\usage{
major(x,y)
}
\arguments{
\item{x,y}{vectors containing non-negative elements (with same length and same mean)}
}
\value{
  logical. TRUE if x majorizes y (x >=(M) y), FALSE if not.
}
\details{even if \code{x} and \code{y} are comparable (i.e. have same length
and same mean) it is possible that neither x majorizes y nor y majorizes x.}
\references{Marshall / Olkin: Inequalities: Theory of Majorization and Its
Applications, New York 1979 (Academic Press)}

\seealso{\code{\link{Lc}}}

\examples{
# generate vectors (of incomes)
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
y <- c(841, 2063, 2445, 3438, 4437, 5401, 6392, 8304, 11304, 21961)
# test whether x majorizes y (TRUE, because y is result of
# Pigou-Dalton-transfers)
major(x,y)
}
\keyword{logic}

