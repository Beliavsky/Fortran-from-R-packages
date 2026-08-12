\name{cl_pam}
\alias{cl_pam}
\title{K-Medoids Partitions of Clusterings}
\description{
  Compute \eqn{k}-medoids partitions of clusterings.
}
\usage{
cl_pam(x, k, method = "euclidean", solver = c("pam", "kmedoids"))
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies, or something
    coercible to that (see \code{\link{cl_ensemble}}).}
  \item{k}{an integer giving the number of classes to be used in the
    partition.}
  \item{method}{a character string or a function, as for argument
    \code{method} of function \code{\link{cl_dissimilarity}}.}
  \item{solver}{a character string indicating the \eqn{k}-medoids solver
    to be employed.  May be abbreviated.  If \code{"pam"} (default), the
    Partitioning Around Medoids
    \bibcitep{|Kaufman+Rousseeuw:1990|chapter 2}
    heuristic \code{\link[cluster]{pam}} of package \pkg{cluster} is
    used.  Otherwise, the exact algorithm of \code{\link{kmedoids}} is
    employed.}
}
\value{
  An object of class \code{"cl_pam"} representing the obtained
  \dQuote{secondary} partition, which is a list with the following
  components.
  \item{cluster}{the class ids of the partition.}
  \item{medoid_ids}{the indices of the medoids.}
  \item{prototypes}{a cluster ensemble with the \eqn{k} prototypes
    (medoids).}
  \item{criterion}{the value of the criterion function of the
    partition.}
  \item{description}{a character string indicating the dissimilarity
    method employed.}
}  
\details{
  An optimal \eqn{k}-medoids partition of the given cluster ensemble is
  defined as a partition of the objects \eqn{x_i} (the elements of the
  ensemble) into \eqn{k} classes \eqn{C_1, \ldots, C_k} such that the
  criterion function
  \eqn{L = \sum_{l=1}^k \min_{j \in C_l} \sum_{i \in C_l} d(x_i, x_j)} 
  is minimized.

  Such secondary partitions \bibcitep{e.g.,|Gordon+Vichi:1998|}
  are obtained by
  computing the dissimilarities \eqn{d} of the objects in the ensemble
  for the given dissimilarity method, and applying a dissimilarity-based
  \eqn{k}-medoids solver to \eqn{d}.
}
\references{
  \bibshow{Kaufman+Rousseeuw:1990, Gordon+Vichi:1998}
}
\seealso{
  \code{\link{cl_pclust}} for more general prototype-based partitions of
  clusterings.
}
\examples{
data("Kinship82")
party <- cl_pam(Kinship82, 3, "symdiff")
## Compare results with tables 5 and 6 in Gordon & Vichi (1998).
party
lapply(cl_prototypes(party), cl_classes)
table(cl_class_ids(party))
}
\keyword{cluster}

