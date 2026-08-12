\name{cl_classes}
\alias{cl_classes}
\title{Cluster Classes}
\description{
  Extract the classes in a partition or hierarchy.
}
\usage{
cl_classes(x)
}
\arguments{
  \item{x}{an R object representing a partition or hierarchy of
    objects.}
}
\value{
  A list inheriting from \code{"cl_classes_of_objects"} of vectors
  indicating the classes.
}
\details{
  For partitions, the classes are the equivalence classes
  (\dQuote{clusters}) of the partition; for soft partitions, the classes
  of the nearest hard partition are used.

  For hierarchies represented by trees, the classes are the sets of
  objects corresponding to (joined at or split by) the nodes of the
  tree.
}
\keyword{cluster}
