\name{lattice}
\encoding{UTF-8}
\alias{cl_meet}
\alias{cl_join}
\alias{Ops.cl_partition}
\alias{Summary.cl_partition}
\alias{Ops.cl_dendrogram}
\alias{Ops.cl_hierarchy}
\alias{Summary.cl_hierarchy}
\title{Cluster Lattices}
\description{
  Computations on the lattice of all (hard) partitions, or the lattice
  of all dendrograms, or the meet semilattice of all hierarchies
  (\eqn{n}-trees) of/on a set of objects: meet, join, and comparisons.
}
\usage{
cl_meet(x, y)
cl_join(x, y)
}
\arguments{
  \item{x}{an ensemble of partitions or dendrograms or hierarchies, or
    an R object representing a partition or dendrogram or hierarchy.}
  \item{y}{an R object representing a partition or dendrogram or
    hierarchy.  Ignored if \code{x} is an ensemble.}
}
\details{
  For a given finite set of objects \eqn{X}, the set \eqn{H(X)} of all
  (hard) partitions of \eqn{X} can be partially ordered by defining a
  partition \eqn{P} to be \dQuote{finer} than a partition \eqn{Q}, i.e.,
  \eqn{P \le Q}, if each class of \eqn{P} is contained in some class of
  \eqn{Q}.  With this partial order, \eqn{H(X)} becomes a bounded
  \dfn{lattice}, with intersection and union of two elements given by
  their greatest lower bound (\dfn{meet}) and their least upper bound
  (\dfn{join}), respectively.

  Specifically, the meet of two partitions computed by \code{cl_meet} is
  the partition obtained by intersecting the classes of the partitions;
  the classes of the join computed by \code{cl_join} are obtained by
  joining all elements in the same class in at least one of the
  partitions.  Obviously, the least and greatest elements of the
  partition lattice are the partitions where each object is in a single
  class (sometimes referred to as the \dQuote{splitter} partition) or in
  the same class (the \dQuote{lumper} partition), respectively.  Meet
  and join of an arbitrary number of partitions can be defined
  recursively.

  In addition to computing the meet and join, the comparison operations
  corresponding to the above partial order as well as \code{min},
  \code{max}, and \code{range} are available at least for R objects
  representing partitions inheriting from \code{"\link{cl_partition}"}.
  The summary methods give the meet and join of the given partitions
  (for \code{min} and \code{max}), or a partition ensemble with the meet
  and join (for \code{range}).

  If the partitions specified by \code{x} and \code{y} are soft
  partitions, the corresponding nearest hard partitions are used.
  Future versions may optionally provide suitable \dQuote{soft} (fuzzy)
  extensions for computing meets and joins.

  The set of all dendrograms on \eqn{X} can be ordered using pointwise
  inequality of the associated ultrametric dissimilarities: i.e., if
  \eqn{D} and \eqn{E} are the dendrograms with ultrametrics \eqn{u} and
  \eqn{v}, respectively, then \eqn{D \le E} if \eqn{u_{ij} \le v_{ij}}
  for all pairs \eqn{(i, j)} of objects.  This again yields a lattice
  (of dendrograms).  The join of \eqn{D} and \eqn{E} is the dendrogram
  with ultrametrics given by \eqn{\max(u_{ij}, v_{ij})} (as this gives
  an ultrametric); the meet is the dendrogram with the maximal
  ultrametric dominated by \eqn{\min(u_{ij}, v_{ij})}, and can be
  obtained by applying single linkage hierarchical clustering to the
  minima.

  The set of all hierarchies on \eqn{X} can be ordered by set-wise
  inclusion of the classes: i.e., if \eqn{H} and \eqn{G} are two
  hierarchies, then \eqn{H \le G} if all classes of \eqn{H} are also
  classes of \eqn{G}.  This yields a meet semilattice, with meet given
  by the classes contained in both hierarchies.  The join only exists if
  the union of the classes is a hierarchy.

  In each case, a modular semilattice is obtained, which allows for a
  natural metrization via least element (semi)lattice move distances,
  see \bibcitet{Barthélemy+Leclerc+Monjardet:1986}.
  These latticial metrics
  are given by the BA/C (partitions), Manhattan (dendrograms), and
  symdiff (hierarchies) dissimilarities, respectively (see
  \code{\link{cl_dissimilarity}}).
}
\value{
  For \code{cl_meet} and \code{cl_join}, an object of class
  \code{"\link{cl_partition}"} or \code{"\link{cl_dendrogram}"} with the
  class ids or ultrametric dissimilarities of the meet and join of the
  partitions or dendrograms, respectively.
}
\references{
  \bibshow{Barthélemy+Leclerc+Monjardet:1986}
}
\examples{
## Two simple partitions of 7 objects.
A <- as.cl_partition(c(1, 1, 2, 3, 3, 5, 5))
B <- as.cl_partition(c(1, 2, 2, 3, 4, 5, 5))
## These disagree on objects 1-3, A splits objects 4 and 5 into
## separate classes.  Objects 6 and 7 are always in the same class.
(A <= B) || (B <= A)
## (Neither partition is finer than the other.)
cl_meet(A, B)
cl_join(A, B)
## Meeting with the lumper (greatest) or joining with the splitter
## (least) partition does not make a difference: 
C_lumper <- as.cl_partition(rep(1, n_of_objects(A)))
cl_meet(cl_ensemble(A, B, C_lumper))
C_splitter <- as.cl_partition(seq_len(n_of_objects(A)))
cl_join(cl_ensemble(A, B, C_splitter))
## Another way of computing the join:
range(A, B, C_splitter)$max
}
\keyword{cluster}
