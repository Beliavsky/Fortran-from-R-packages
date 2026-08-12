\name{cl_ultrametric}
\alias{cl_ultrametric}
\alias{as.cl_ultrametric}
\title{Ultrametrics of Hierarchies}
\description{
  Compute the ultrametric distances for objects representing (total
  indexed) hierarchies.
}
\usage{
cl_ultrametric(x, size = NULL, labels = NULL)
as.cl_ultrametric(x)
}
\arguments{
  \item{x}{an R object representing a (total indexed) hierarchy of
    objects.}
  \item{size}{an integer giving the number of objects in the hierarchy.}
  \item{labels}{a character vector giving the names of the objects in
    the hierarchy.}
}
\value{
  An object of class \code{"cl_ultrametric"} containing the ultrametric
  distances.
}  
\details{
  If \code{x} is not an ultrametric or a hierarchy with an ultrametric
  representation, \code{cl_ultrametric} uses
  \code{\link[stats]{cophenetic}} to obtain the ultrametric (also known
  as cophenetic) distances from the hierarchy, which in turn by default
  calls the S3 generic \code{\link[stats]{as.hclust}} on the hierarchy.
  Support for a class which represents hierarchies can thus be added by 
  providing \code{as.hclust} methods for this class.  In R 2.1.0 or
  better, \code{cophenetic} is an S3 generic as well, and one can also
  more directly provide methods for this if necessary.

  \code{as.cl_ultrametric} is a generic function which can be used for
  coercing \emph{raw} (non-classed) ultrametrics, represented as numeric
  vectors (of the lower-half entries) or numeric matrices, to
  ultrametric objects.

  Ultrametric objects are implemented as symmetric proximity objects
  with a dissimilarity interpretation so that self-proximities are zero,
  and inherit from classes \code{"\link{cl_dissimilarity}"} and
  \code{"cl_proximity"}.  See section \bold{Details} in the
  documentation for \code{\link{cl_dissimilarity}} for implications.

  Ultrametric objects can also be coerced to classes
  \code{"\link[stats]{dendrogram}"} and
  \code{"\link[stats]{hclust}"}, and hence in particular use the
  \code{plot} methods for these classes.  By default, plotting an
  ultrametric object uses the plot method for dendrograms.
}
\seealso{
  \code{\link{is.cl_hierarchy}}
}
\examples{
hc <- hclust(dist(USArrests))
u <- cl_ultrametric(hc)
## Subscripting.
u[1 : 5, 1 : 5]
u[1 : 5, 6 : 7]
## Plotting.
plot(u)
}
\keyword{cluster}
