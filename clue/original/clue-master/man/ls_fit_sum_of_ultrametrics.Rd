\name{ls_fit_sum_of_ultrametrics}
\alias{ls_fit_sum_of_ultrametrics}
\title{Least Squares Fit of Sums of Ultrametrics to Dissimilarities}
\description{
  Find a sequence of ultrametrics with sum minimizing square distance
  (Euclidean dissimilarity) to a given dissimilarity object.
}
\usage{
ls_fit_sum_of_ultrametrics(x, nterms = 1, weights = 1,
                           control = list())
}
\arguments{
  \item{x}{a dissimilarity object inheriting from or coercible to class
    \code{"\link{dist}"}.}
  \item{nterms}{an integer giving the number of ultrametrics to be
    fitted.}
  \item{weights}{a numeric vector or matrix with non-negative weights
    for obtaining a weighted least squares fit.  If a matrix, its
    numbers of rows and columns must be the same as the number of
    objects in \code{x}, and the lower diagonal part is used.
    Otherwise, it is recycled to the number of elements in \code{x}.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\details{
  The problem to be solved is minimizing the criterion function
  \deqn{L(u(1), \dots, u(n)) =
    \sum_{i,j} w_{ij} \left(x_{ij} - \sum_{k=1}^n u_{ij}(k)\right)^2}{
    L(u(1), \dots, u(n)) =
    \sum_{i,j} w_{ij} \left(x_{ij} - \sum_{k=1}^n u_{ij}(k)\right)^2}
  over all \eqn{u(1), \ldots, u(n)} satisfying the ultrametric
  constraints.

  We provide an implementation of the iterative heuristic suggested in
  \bibcitet{Carroll+Pruzansky:1980} which in each step \eqn{t} sequentially
  refits the \eqn{u(k)} as the least squares ultrametric fit to the
  \dQuote{residuals} \eqn{x - \sum_{l \ne k} u(l)} using
  \code{\link{ls_fit_ultrametric}}.

  Available control parameters include
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      iteration steps to be performed.
      Defaults to 100.}
    \item{\code{eps}}{a nonnegative number controlling the iteration,
      which stops when the maximal change in all \eqn{u(k)} is less than
      \code{eps}.
      Defaults to \eqn{10^{-6}}.}
    \item{\code{reltol}}{the relative convergence tolerance.  Iteration
      stops when the relative change in the criterion function is less
      than \code{reltol}.
      Defaults to \eqn{10^{-6}}.}
    \item{\code{method}}{a character string indicating the fitting
      method to be employed by the individual least squares fits.}
    \item{\code{control}}{a list of control parameters to be used by the
      method of \code{\link{ls_fit_ultrametric}} employed.  By default,
      if the \acronym{SUMT} method method is used, 10 inner
      \acronym{SUMT} runs are performed for each refitting.}
  }

  It should be noted that the method used is a heuristic which can not
  be guaranteed to find the global minimum.
}
\value{
  A list of objects of class \code{"\link{cl_ultrametric}"} containing
  the fitted ultrametric distances.
}
\references{
  \bibshow{Carroll+Pruzansky:1980}
}
\keyword{cluster}
\keyword{optimize}
