\name{volume.ellipsoid}
\alias{volume}
\alias{volume.ellipsoid}
\title{Compute the Volume (of an Ellipsoid)}
\description{
  Compute the volume of geometric \R object.
  This is a generic function and has a method for \code{ellipsoid} objects
  (typically resulting from \code{\link{ellipsoidhull}()}.
}
\usage{
volume(object, \dots)
\method{volume}{ellipsoid}(object, log = FALSE, \dots)
}
\arguments{
  \item{object}{an \R object the volume of which is wanted; for the
    \code{ellipsoid} method, an object of that class (see
    \code{\link{ellipsoidhull}} or the example below).}
  \item{log}{\code{\link{logical}} indicating if the volume should be
    returned in log scale.  Maybe needed in largish dimensions.}
  \item{\dots}{potential further arguments of methods, e.g. \code{log}.}
}
\value{
  a number, the volume \eqn{V} (or \eqn{\log(V)} if \code{log = TRUE}) of
  the given \code{object}.
}
\author{
  Martin Maechler (2002, extracting from former \code{\link{clusplot}} code);
  Keefe Murphy (2019) provided code for dimensions \eqn{d > 2}.
}
\seealso{\code{\link{ellipsoidhull}} for spanning ellipsoid computation.}
\examples{
## example(ellipsoidhull) # which defines 'ellipsoid' object <namefoo>

myEl <- structure(list(cov = rbind(c(3,1),1:2), loc = c(0,0), d2 = 10),
                   class = "ellipsoid")
volume(myEl)# i.e. "area" here (d = 2)
myEl # also mentions the "volume"

set.seed(1)
d5 <- matrix(rt(500, df=3), 100,5)
e5 <- ellipsoidhull(d5)
}
\keyword{utilities}
