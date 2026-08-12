\name{partition}
\alias{cl_partition}		% class ...
\alias{is.cl_partition}
\alias{is.cl_hard_partition}
\alias{is.cl_soft_partition}
\alias{cl_hard_partition}	% class ...
\alias{as.cl_partition}
\alias{as.cl_hard_partition}
\title{Partitions}
\description{
  Determine whether an R object represents a partition of objects, or
  coerce to an R object representing such.}
\usage{
is.cl_partition(x)
is.cl_hard_partition(x)
is.cl_soft_partition(x)

as.cl_partition(x)
as.cl_hard_partition(x)
}
\arguments{
  \item{x}{an R object.}
}
\value{
  For the testing functions, a logical indicating whether the given
  object represents a clustering of objects of the respective kind.

  For the coercion functions, a container object inheriting from
  \code{"cl_partition"}, with a suitable representation of the partition
  given by \code{x}.
}
\details{
  \code{is.cl_partition} and \code{is.cl_hard_partition} are generic
  functions.

  The methods provided in package \pkg{clue} handle the partitions
  obtained from clustering functions in the base R distribution, as well
  as packages \pkg{RWeka}, \pkg{cba}, \pkg{cclust}, \pkg{cluster},
  \pkg{e1071}, \pkg{flexclust}, \pkg{flexmix}, \pkg{kernlab},
  \pkg{mclust}, \pkg{movMF} and \pkg{skmeans} (and of course, \pkg{clue}
  itself).

  \code{is.cl_soft_partition} gives true iff \code{is.cl_partition} is
  true and \code{is.cl_hard_partition} is false.

  \code{as.cl_partition} returns an object of class
  \code{"cl_partition"} \dQuote{containing} the given object \code{x} if
  this already represents a partition (i.e., \code{is.cl_partition(x)}
  is true), or the memberships obtained from \code{x} via
  \code{\link{as.cl_membership}}.

  \code{as.cl_hard_partition(x)} returns an object which has class
  \code{"cl_hard_partition"} and inherits from \code{"cl_partition"},
  and contains \code{x} if it already represents a hard partition (i.e.,
  provided that \code{is.cl_hard_partition(x)} is true), or the class
  ids obtained from \code{x}, using \code{x} if this is an atomic vector
  of raw class ids, or, if \code{x} represents a soft partition or is a
  raw matrix of membership values, using the class ids of the
  \emph{nearest hard partition}, defined by taking the class ids of the
  (first) maximal membership values.

  Conceptually, partitions and hard partitions are \emph{virtual}
  classes, allowing for a variety of representations.

  There are group methods for comparing partitions and computing their
  minimum, maximum, and range based on the meet and join operations, see
  \code{\link{cl_meet}}.
}
\examples{
data("Cassini")
pcl <- kmeans(Cassini$x, 3)
is.cl_partition(pcl)
is.cl_hard_partition(pcl)
is.cl_soft_partition(pcl)
}
\keyword{cluster}
