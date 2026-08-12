\name{l1_fit_ultrametric}
\alias{l1_fit_ultrametric}
\title{Least Absolute Deviation Fit of Ultrametrics to Dissimilarities}
\description{
  Find the ultrametric with minimal absolute distance (Manhattan
  dissimilarity) to a given dissimilarity object.
}
\usage{
l1_fit_ultrametric(x, method = c("SUMT", "IRIP"), weights = 1,
                   control = list())
}
\arguments{
  \item{x}{a dissimilarity object inheriting from or coercible to class
    \code{"\link{dist}"}.}
  \item{method}{a character string indicating the fitting method to be
    employed.  Must be one of \code{"SUMT"} (default) or \code{"IRIP"},
    or a unique abbreviation thereof.}
  \item{weights}{a numeric vector or matrix with non-negative weights
    for obtaining a weighted least squares fit.  If a matrix, its
    numbers of rows and columns must be the same as the number of
    objects in \code{x}, and the lower diagonal part is used.
    Otherwise, it is recycled to the number of elements in \code{x}.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  An object of class \code{"\link{cl_ultrametric}"} containing the
  fitted ultrametric distances.
}
\details{
  The problem to be solved is minimizing
  \deqn{L(u) = \sum_{i,j} w_{ij} |x_{ij} - u_{ij}|}
  over all \eqn{u} satisfying the ultrametric constraints (i.e., for all
  \eqn{i, j, k}, \eqn{u_{ij} \le \max(u_{ik}, u_{jk})}).  This problem
  is known to be NP hard
  \bibcitep{Křivánek+Morávek:1986}.

  We provide two heuristics for solving this problem.

  Method \code{"SUMT"} implements a \acronym{SUMT} (Sequential
  Unconstrained Minimization Technique, see \code{\link{sumt}}) approach
  using the sign function for the gradients of the absolute value
  function.

  Available control parameters are \code{method}, \code{control},
  \code{eps}, \code{q}, and \code{verbose}, which have the same roles as
  for \code{\link{sumt}}, and the following.

  \describe{
    \item{\code{nruns}}{an integer giving the number of runs to be
      performed.  Defaults to 1.}
    \item{\code{start}}{a single dissimilarity, or a list of
      dissimilarities to be employed as starting values.}
  }

  Method \code{"IRIP"} implements a variant of the Iteratively
  Reweighted Iterative Projection approach of \bibcitet{Smith:2001},
  which
  attempts to solve the \eqn{L_1} problem via a sequence of weighted
  \eqn{L_2} problems, determining \eqn{u(t+1)} by minimizing the
  criterion function
  \deqn{\sum_{i,j} w_{ij}
    (x_{ij} - u_{ij})^2 / \max(|x_{ij} - u_{ij}(t)|, m)}
  with \eqn{m} a \dQuote{small} non-zero value to avoid zero divisors.
  We use the \acronym{SUMT} method of \code{\link{ls_fit_ultrametric}}
  for solving the weighted least squares problems.

  Available control parameters are as follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      iteration steps to be performed.
      Defaults to 100.}
    \item{\code{eps}}{a nonnegative number controlling the iteration,
      which stops when the maximal change in \eqn{u} is less than
      \code{eps}.
      Defaults to \eqn{10^{-6}}.}
    \item{\code{reltol}}{the relative convergence tolerance.  Iteration
      stops when the relative change in the \eqn{L_1} criterion is less
      than \code{reltol}.
      Defaults to \eqn{10^{-6}}.}
    \item{\code{MIN}}{the cutoff \eqn{m}.  Defaults to \eqn{10^{-3}}.}
    \item{\code{start}}{a dissimilarity object to be used as the
      starting value for \eqn{u}.}
    \item{\code{control}}{a list of control parameters to be used by the
      method of \code{\link{ls_fit_ultrametric}} employed for solving
      the weighted \eqn{L_2} problems.}
  }
    
  One may need to adjust the default control parameters to achieve
  convergence.

  It should be noted that all methods are heuristics which can not be
  guaranteed to find the global minimum.  
}
\seealso{
  \code{\link{cl_consensus}} for computing least absolute deviation
  (Manhattan) consensus hierarchies;
  \code{\link{ls_fit_ultrametric}}.
}
\references{
  \bibshow{Křivánek+Morávek:1986, Smith:2001}
}  
\keyword{cluster}
\keyword{optimize}
