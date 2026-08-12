\name{fit_ultrametric_target}
\alias{ls_fit_ultrametric_target}
\alias{l1_fit_ultrametric_target}
\title{Fit Dissimilarities to a Hierarchy}
\description{
  Find the ultrametric from a target equivalence class of hierarchies
  which minimizes weighted Euclidean or Manhattan dissimilarity to a
  given dissimilarity object.
}
\usage{
ls_fit_ultrametric_target(x, y, weights = 1)
l1_fit_ultrametric_target(x, y, weights = 1)
}
\arguments{
  \item{x}{a dissimilarity object inheriting from class
    \code{"\link{dist}"}.}
  \item{y}{a target hierarchy.}
  \item{weights}{a numeric vector or matrix with non-negative weights
    for obtaining a weighted fit.  If a matrix, its numbers of rows and
    columns must be the same as the number of objects in \code{x}.
    Otherwise, it is recycled to the number of elements in \code{x}.}
}
\value{
  An object of class \code{"\link{cl_ultrametric}"} containing the
  optimal ultrametric distances.
}
\details{
  The target equivalence class consists of all dendrograms for which the
  corresponding \eqn{n}-trees are the same as the one corresponding to
  \code{y}.  I.e., all splits are the same as for \code{y}, and
  optimization is over the height of the splits.

  The criterion function to be optimized over all ultrametrics from the
  equivalence class is \eqn{\sum w_{ij} |x_{ij} - u_{ij}|^p}, where
  \eqn{p = 2} in the Euclidean and \eqn{p = 1} in the Manhattan case,
  respectively.

  The optimum can be computed as follows.  Suppose split \eqn{s} joins
  object classes \eqn{A} and \eqn{B}.  As the ultrametric
  dissimilarities of all objects in \eqn{A} to all objects in \eqn{B}
  must be the same value, say, \eqn{u_{A,B} = u_s}, the contribution
  from the split to the criterion function is of the form
  \eqn{f_s(u_s) = \sum_{i \in A, j \in B} w_{ij} |x_{ij} - u_s|^p}.
  We need to minimize \eqn{\sum_s f_s(u_s)} under the constraint that
  the \eqn{u_s} form a non-decreasing sequence, which is accomplished by
  using the Pool Adjacent Violator Algorithm (\acronym{PAVA}) using the
  weighted mean (\eqn{p = 2}) or weighted median (\eqn{p = 1}) for
  solving the blockwise optimization problems.
}
\seealso{
  \code{\link{ls_fit_ultrametric}} for finding the ultrametric
  minimizing Euclidean dissimilarity (without fixing the splits).
}
\examples{
data("Phonemes")
## Note that the Phonemes data set has the consonant misclassification
## probabilities, i.e., the similarities between the phonemes.
d <- as.dist(1 - Phonemes)
## Find the maximal dominated and miminal dominating ultrametrics by
## hclust() with single and complete linkage:
y1 <- hclust(d, "single")
y2 <- hclust(d, "complete")
## Note that these are quite different:
cl_dissimilarity(y1, y2, "gamma")
## Now find the L2 optimal members of the respective dendrogram
## equivalence classes.
u1 <- ls_fit_ultrametric_target(d, y1)
u2 <- ls_fit_ultrametric_target(d, y2)
## Compute the L2 optimal ultrametric approximation to d.
u <- ls_fit_ultrametric(d)
## And compare ...
cl_dissimilarity(cl_ensemble(Opt = u, Single = u1, Complete = u2), d)
## The solution obtained via complete linkage is quite close:
cl_agreement(u2, u, "cophenetic")
}
\keyword{cluster}
\keyword{optimize}
