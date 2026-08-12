\name{kmedoids}
\alias{kmedoids}
\title{K-Medoids Clustering}
\description{
  Compute a \eqn{k}-medoids partition of a dissimilarity object.
}
\usage{
kmedoids(x, k)
}
\arguments{
  \item{x}{a dissimilarity object inheriting from class
    \code{"\link{dist}"}, or a square matrix of pairwise
    object-to-object dissimilarity values.}
  \item{k}{an integer giving the number of classes to be used in the
    partition.}
}
\value{
  An object of class \code{"kmedoids"} representing the obtained
  partition, which is a list with the following components.
  \item{cluster}{the class ids of the partition.}
  \item{medoid_ids}{the indices of the medoids.}
  \item{criterion}{the value of the criterion function of the
    partition.}
}
\details{
  Let \eqn{d} denote the pairwise object-to-object dissimilarity matrix
  corresponding to \code{x}.  A \eqn{k}-medoids partition of \code{x} is
  defined as a partition of the numbers from 1 to \eqn{n}, the number of
  objects in \code{x}, into \eqn{k} classes \eqn{C_1, \ldots, C_k} such
  that the criterion function
  \eqn{L = \sum_l \min_{j \in C_l} \sum_{i \in C_l} d_{ij}}
  is minimized.

  This is an NP-hard optimization problem.  PAM (Partitioning Around
  Medoids, see \bibcitet{|Kaufman+Rousseeuw:1990|chapter 2})
  is a very popular
  heuristic for obtaining optimal \eqn{k}-medoids partitions, and
  provided by \code{\link[cluster]{pam}} in package \pkg{cluster}.

  \code{kmedoids} is an exact algorithm based on a binary linear
  programming formulation of the optimization problem
  \bibcitep{e.g.,|Gordon+Vichi:1998|[P4']},
  using \code{\link[lpSolve]{lp}} from package
  \pkg{lpSolve} as solver.  Depending on available hardware resources
  (the number of constraints of the program is of the order \eqn{n^2}),
  it may not be possible to obtain a solution.
}
\references{
  \bibshow{Kaufman+Rousseeuw:1990, Gordon+Vichi:1998}
}  
\keyword{cluster}
\keyword{optimize}
