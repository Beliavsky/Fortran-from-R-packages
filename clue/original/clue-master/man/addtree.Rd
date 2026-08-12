\name{addtree}
\encoding{UTF-8}
\alias{as.cl_addtree}
\title{Additive Tree Distances}
\description{
  Objects representing additive tree distances.
}
\usage{
as.cl_addtree(x)
}
\arguments{
  \item{x}{an R object representing additive tree distances.}
}
\value{
  An object of class \code{"cl_addtree"} containing the additive
  tree distances.
}
\details{
  Additive tree distances are object dissimilarities \eqn{d} satisfying
  the so-called \emph{additive tree conditions}, also known as
  \emph{four-point conditions} \eqn{d_{ij} + d_{kl} \le \max(d_{ik} +
    d_{jl}, d_{il} + d_{jk})} for all quadruples \eqn{i, j, k, l}.
  Equivalently, for each such quadruple, the largest two values of the
  sums \eqn{d_{ij} + d_{kl}}, \eqn{d_{ik} + d_{jl}}, and \eqn{d_{il} +
    d_{jk}} must be equal.
  Centroid distances are additive tree distances where the inequalities
  in the four-point conditions are strengthened to equalities (such that
  all three sums are equal), and can be represented as \eqn{d_{ij} = g_i
    + g_j}, i.e., as sums of distances from a \dQuote{centroid}.
  See, e.g., \bibcitet{Barthélemy+Guénoche:1991} for more details on
  additive tree distances.

  \code{as.cl_addtree} is a generic function.  Its default method can
  handle objects representing ultrametric distances and raw additive
  distance matrices.  In addition, there is a method for coercing
  objects of class \code{"\link[ape:as.phylo]{phylo}"} from package
  \pkg{ape}.

  Functions \code{\link{ls_fit_addtree}} and
  \code{\link{ls_fit_centroid}} can be used to find the additive tree
  distance or centroid distance minimizing least squares distance
  (Euclidean dissimilarity) to a given dissimilarity object.

  There is a \code{\link{plot}} method for additive tree distances.
}
\references{
  \bibshow{Barthélemy+Guénoche:1991}
}
\keyword{cluster}
