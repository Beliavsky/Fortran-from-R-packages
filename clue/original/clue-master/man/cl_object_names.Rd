\name{cl_object_names}
\alias{cl_object_names}
\title{Find Object Names}
\description{
  Find the names of the objects from which a taxonomy (partition or
  hierarchy) or proximity was obtained.
}
\usage{
cl_object_names(x)
}
\arguments{
  \item{x}{an \R object representing a taxonomy or proximity.}
}
\value{
  A character vector of length \code{\link{n_of_objects}(x)} in case the
  names of the objects could be determined, or \code{NULL}.
}
\details{
  This is a generic function.

  The methods provided in package \pkg{clue} handle the partitions and
  hierarchies obtained from clustering functions in the base R
  distribution, as well as packages \pkg{RWeka}, \pkg{ape}, \pkg{cba},
  \pkg{cclust}, \pkg{cluster}, \pkg{e1071}, \pkg{flexclust},
  \pkg{flexmix}, \pkg{kernlab}, \pkg{mclust}, \pkg{movMF} and
  \pkg{skmeans} (and of course, \pkg{clue} itself), in as much as
  possible.

  There is also a method for object dissimilarities which inherit from
  class \code{"\link{dist}"}.
}
\keyword{cluster}
