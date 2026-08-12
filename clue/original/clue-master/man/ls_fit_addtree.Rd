\name{ls_fit_addtree}
\encoding{UTF-8}
\alias{ls_fit_addtree}
\alias{ls_fit_centroid}
\title{Least Squares Fit of Additive Tree Distances to Dissimilarities}
\description{
  Find the additive tree distance or centroid distance minimizing least
  squares distance (Euclidean dissimilarity) to a given dissimilarity
  object.
}
\usage{
ls_fit_addtree(x, method = c("SUMT", "IP", "IR"), weights = 1,
               control = list())
ls_fit_centroid(x)
}
\arguments{
  \item{x}{a dissimilarity object inheriting from class
    \code{"\link{dist}"}.}
  \item{method}{a character string indicating the fitting method to be
    employed.  Must be one of \code{"SUMT"} (default), \code{"IP"}, or
    \code{"IR"}, or a unique abbreviation thereof.}
  \item{weights}{a numeric vector or matrix with non-negative weights
    for obtaining a weighted least squares fit.  If a matrix, its
    numbers of rows and columns must be the same as the number of
    objects in \code{x}, and the lower diagonal part is used.
    Otherwise, it is recycled to the number of elements in \code{x}.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  An object of class \code{"cl_addtree"} containing the optimal additive
  tree distances.
}
\details{
  See \code{\link{as.cl_addtree}} for details on additive tree distances
  and centroid distances.
  
  With \eqn{L(d) = \sum w_{ij} (x_{ij} - d_{ij})^2}, the problem to be
  solved by \code{ls_fit_addtree} is minimizing \eqn{L} over all
  additive tree distances \eqn{d}.  This problem is known to be NP
  hard.

  We provide three heuristics for solving this problem.

  Method \code{"SUMT"} implements the \acronym{SUMT}
  \bibcitep{Sequential Unconstrained Minimization Technique|Fiacco+McCormick:1968|}
  approach of \bibcitet{De_Soete:1983}.
  Incomplete dissimilarities are currently not supported.

  Methods \code{"IP"} and \code{"IR"} implement the Iterative
  Projection and Iterative Reduction approaches of
  \bibcitet{Hubert+Arabie:1995} and \bibcitet{Roux:1988}, respectively.
  Non-identical weights and
  incomplete dissimilarities are currently not supported.

  See \code{\link{ls_fit_ultrametric}} for details on these methods and
  available control parameters.

  It should be noted that all methods are heuristics which can not be
  guaranteed to find the global minimum.  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.

  \code{ls_fit_centroid} finds the centroid distance \eqn{d} minimizing
  \eqn{L(d)} (currently, only for the case of identical weights).  This
  optimization problem has a closed-form solution.
}
\references{
  \bibshow{Fiacco+McCormick:1968,
    Hubert+Arabie:1995,
    Roux:1988,
    De_Soete:1983}
}
\keyword{cluster}
\keyword{optimize}
