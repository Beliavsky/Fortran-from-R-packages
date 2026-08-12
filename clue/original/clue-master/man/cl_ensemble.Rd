\name{cl_ensemble}
\alias{cl_ensemble}
\alias{as.cl_ensemble}
\alias{is.cl_ensemble}
\title{Cluster Ensembles}
\description{Creation and manipulation of cluster ensembles.}
\usage{
cl_ensemble(..., list = NULL)
as.cl_ensemble(x)
is.cl_ensemble(x)
}
\arguments{
  \item{\dots}{R objects representing clusterings of or dissimilarities
    between the same objects.}
  \item{list}{a list of R objects as in \code{\dots}.}
  \item{x}{for \code{as.cl_ensemble}, an R object as in \code{\dots};
    for \code{is.cl_ensemble}, an arbitrary R object.}
}
\details{
  \code{cl_ensemble} creates \dQuote{cluster ensembles}, which are
  realized as lists of clusterings (or dissimilarities) with additional
  class information, always inheriting from \code{"cl_ensemble"}.  All
  elements of the ensemble must have the same number of objects.

  If all elements are partitions, the ensemble has class
  \code{"cl_partition_ensemble"};
  if all elements are dendrograms, it has class
  \code{"cl_dendrogram_ensemble"} and inherits from
  \code{"cl_hierarchy_ensemble"};
  if all elements are hierarchies (but not always dendrograms), it has
  class \code{"cl_hierarchy_ensemble"}.
  Note that empty or \dQuote{mixed} ensembles cannot be categorized
  according to the kind of elements they contain, and hence only have
  class \code{"cl_ensemble"}.

  The list representation makes it possible to use \code{lapply} for
  computations on the individual clusterings in (i.e., the components
  of) a cluster ensemble.

  Available methods for cluster ensembles include those for
  subscripting, \code{c}, \code{rep}, and \code{print}.  There is also a
  \code{plot} method for ensembles for which all elements can be plotted
  (currently, additive trees, dendrograms and ultrametrics).
}
\value{
  \code{cl_ensemble} returns a list of the given clusterings or
  dissimilarities, with additional class information (see
  \bold{Details}).
}
\examples{
d <- dist(USArrests)
hclust_methods <-
    c("ward", "single", "complete", "average", "mcquitty")
hclust_results <- lapply(hclust_methods, function(m) hclust(d, m))
names(hclust_results) <- hclust_methods 
## Now create an ensemble from the results.
hens <- cl_ensemble(list = hclust_results)
hens
## Subscripting.
hens[1 : 3]
## Replication.
rep(hens, 3)
## Plotting.
plot(hens, main = names(hens))
## And continue to analyze the ensemble, e.g.
round(cl_dissimilarity(hens, method = "gamma"), 4)
}
\keyword{cluster}
