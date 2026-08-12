\name{cl_prototypes}
\alias{cl_prototypes}
\title{Partition Prototypes}
\description{
  Determine prototypes for the classes of an R object representing a
  partition.
}
\usage{
cl_prototypes(x)
}
\arguments{
  \item{x}{an R object representing a partition of objects.}
}
\details{
  Many partitioning methods are based on prototypes (\dQuote{centers},
  \dQuote{centroids}, \dQuote{medoids}, \dots).  In typical cases, these
  are points in the feature space for the measurements on the objects to
  be partitioned, such that one can quantify the distance between the
  objects and the prototypes, and, e.g., classify objects to their
  closest prototype.

  This is a generic function.  The methods provided in package
  \pkg{clue} handle the partitions obtained from clustering functions in
  the base R distribution, as well as packages \pkg{cba}, \pkg{cclust},
  \pkg{cluster}, \pkg{e1071}, \pkg{flexclust}, \pkg{kernlab}, and
  \pkg{mclust} (and of course, \pkg{clue} itself).
}
\examples{
## Show how prototypes ("centers") vary across k-means runs on
## bootstrap samples from the Cassini data.
data("Cassini")
nr <- NROW(Cassini$x)
out <- replicate(50,
                 { kmeans(Cassini$x[sample(nr, replace = TRUE), ], 3) },
                 simplify = FALSE)
## Plot the data points in light gray, and the prototypes found.
plot(Cassini$x, col = gray(0.8))
points(do.call("rbind", lapply(out, cl_prototypes)), pch = 19)
}
\keyword{cluster}
