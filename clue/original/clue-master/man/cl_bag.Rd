\name{cl_bag}
\alias{cl_bag}
\title{Bagging for Clustering}
\description{
  Construct partitions of objects by running a base clustering algorithm
  on bootstrap samples from a given data set, and \dQuote{suitably}
  aggregating these primary partitions.
}
\usage{
cl_bag(x, B, k = NULL, algorithm = "kmeans", parameters = NULL, 
       method = "DFBC1", control = NULL)
}
\arguments{
  \item{x}{the data set of objects to be clustered, as appropriate for
    the base clustering algorithm.}
  \item{B}{an integer giving the number of bootstrap replicates.}
  \item{k}{\code{NULL} (default), or an integer giving the number of
    classes to be used for a partitioning base algorithm.}
  \item{algorithm}{a character string or function specifying the base
    clustering algorithm.}
  \item{parameters}{a named list of additional arguments to be passed to
    the base algorithm.}
  \item{method}{a character string indicating the bagging method to
    use.  Currently, only method \code{"DFBC1"} is available, which
    implements algorithm \emph{BagClust1} in
    \bibcitet{Dudoit+Fridlyand:2003}.}
  \item{control}{a list of control parameters for the aggregation.
    Currently, not used.}
}
\value{
  An R object representing a partition of the objects given in \code{x}.
}
\details{
  Bagging for clustering is really a rather general conceptual framework
  than a specific algorithm.  If the primary partitions generated in the
  bootstrap stage form a cluster ensemble (so that class memberships of
  the objects in \code{x} can be obtained), consensus methods for
  cluster ensembles (as implemented, e.g., in \code{\link{cl_consensus}}
  and \code{\link{cl_medoid}}) can be employed for the aggregation
  stage.  In particular, (possibly new) bagging algorithms can easily be
  realized by directly running \code{\link{cl_consensus}} on the results
  of \code{\link{cl_boot}}.

  In BagClust1, aggregation proceeds by generating a reference partition
  by running the base clustering algorithm on the whole given data set,
  and averaging the ensemble memberships after optimally matching them
  to the reference partition (in fact, by minimizing Euclidean
  dissimilarity, see \code{\link{cl_dissimilarity}}).

  If the base clustering algorithm yields prototypes, aggregation can be
  based on clustering these.  This is the idea underlying the
  \dQuote{Bagged Clustering} algorithm introduced in
  \bibcitet{Leisch:1999} and implemented by function
  \code{\link[e1071]{bclust}} in package \pkg{e1071}.
}
\references{
  \bibshow{Dudoit+Fridlyand:2003, Leisch:1999}
}
\examples{
set.seed(1234)
## Run BagClust1 on the Cassini data.
data("Cassini")
party <- cl_bag(Cassini$x, 50, 3)
plot(Cassini$x, col = cl_class_ids(party), xlab = "", ylab = "")
## Actually, using fuzzy c-means as a base learner works much better:
if(require("e1071", quietly = TRUE)) {
    party <- cl_bag(Cassini$x, 20, 3, algorithm = "cmeans")
    plot(Cassini$x, col = cl_class_ids(party), xlab = "", ylab = "")
}
}
\keyword{cluster}
