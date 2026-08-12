\name{n_of_classes}
\alias{n_of_classes}
\alias{cl_class_ids}
\alias{as.cl_class_ids}
\title{Classes in a Partition}
\description{Determine the number of classes and the class ids in a
  partition of objects.}
\usage{
n_of_classes(x)
cl_class_ids(x)
as.cl_class_ids(x)
}
\arguments{
  \item{x}{an object representing a (hard or soft) partition (for
    \code{n_of_classes} and \code{cl_class_ids}), or raw class ids (for
    \code{as.cl_class_ids}).}
}
\value{
  For \code{n_of_classes}, an integer giving the number of classes in
  the partition.
  
  For \code{cl_class_ids}, a vector of integers with the corresponding
  class ids.  For soft partitions, the class ids returned are those of
  the \emph{nearest hard partition} obtained by taking the class ids of
  the (first) maximal membership values.
}
\details{
  These function are generic functions.

  The methods provided in package \pkg{clue} handle the partitions
  obtained from clustering functions in the base R distribution, as well
  as packages \pkg{RWeka}, \pkg{cba}, \pkg{cclust}, \pkg{cluster},
  \pkg{e1071}, \pkg{flexclust}, \pkg{flexmix}, \pkg{kernlab},
  \pkg{mclust}, \pkg{movMF} and \pkg{skmeans} (and of course, \pkg{clue}
  itself).

  Note that the number of classes is taken as the number of distinct
  class ids actually used in the partition, and may differ from the
  number of columns in a membership matrix representing the partition.

  \code{as.cl_class_ids} can be used for coercing \dQuote{raw} class
  ids (given as atomic vectors) to class id objects.
}
\seealso{
  \code{\link{is.cl_partition}}
}
\examples{
data("Cassini")
party <- kmeans(Cassini$x, 3)
n_of_classes(party)
## A simple confusion matrix:
table(cl_class_ids(party), Cassini$classes)
## For an "oversize" membership matrix representation:
n_of_classes(cl_membership(party, 6))
}
\keyword{cluster}
