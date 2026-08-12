\name{ls_fit_ultrametric}
\encoding{UTF-8}
\alias{ls_fit_ultrametric}
\title{Least Squares Fit of Ultrametrics to Dissimilarities}
\description{
  Find the ultrametric with minimal square distance (Euclidean
  dissimilarity) to given dissimilarity objects.
}
\usage{
ls_fit_ultrametric(x, method = c("SUMT", "IP", "IR"), weights = 1,
                   control = list())
}
\arguments{
  \item{x}{a dissimilarity object inheriting from or coercible to class
    \code{"\link{dist}"}, or an ensemble of such objects.} 
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
  An object of class \code{"\link{cl_ultrametric}"} containing the
  fitted ultrametric distances.
}
\details{
  For a single dissimilarity object \code{x}, the problem to be solved
  is minimizing
  \deqn{L(u) = \sum_{i,j} w_{ij} (x_{ij} - u_{ij})^2}
  over all \eqn{u} satisfying the ultrametric constraints (i.e., for all
  \eqn{i, j, k}, \eqn{u_{ij} \le \max(u_{ik}, u_{jk})}).  This problem
  is known to be NP hard \bibcitep{Křivánek+Morávek:1986}.

  For an ensemble of dissimilarity objects, the criterion function is
  \deqn{L(u) = \sum_b w_b \sum_{i,j} w_{ij} (x_{ij}(b) - u_{ij})^2,}
  where \eqn{w_b} is the weight given to element \eqn{x_b} of the
  ensemble and can be specified via control parameter \code{weights}
  (default: all ones).  This problem reduces to the above basic problem
  with \eqn{x} as the \eqn{w_b}-weighted mean of the \eqn{x_b}.

  We provide three heuristics for solving the basic problem.

  Method \code{"SUMT"} implements the \acronym{SUMT}
  \bibcitep{Sequential Unconstrained Minimization Technique|Fiacco+McCormick:1968|}
  approach of \bibcitet{De_Soete:1984a} which in turn simplifies the suggestions
  in \bibcitet{Carroll+Pruzansky:1980}.  (See \code{\link{sumt}} for more
  information on the \acronym{SUMT} approach.)  We then use a final
  single linkage hierarchical clustering step to ensure that the
  returned object exactly satisfies the ultrametric constraints.  The
  starting value \eqn{u_0} is obtained by \dQuote{random shaking} of the
  given dissimilarity object (if not given).  If there are missing
  values in \code{x}, i.e., the given dissimilarities are
  \emph{incomplete}, we follow a suggestion of \bibcitet{De_Soete:1984a}, imputing
  the missing values by the weighted mean of the non-missing ones, and
  setting the corresponding weights to zero.

  Available control parameters are \code{method}, \code{control},
  \code{eps}, \code{q}, and \code{verbose}, which have the same roles as
  for \code{\link{sumt}}, and the following.
  
  \describe{
    \item{\code{nruns}}{an integer giving the number of runs to be
      performed.  Defaults to 1.}
    \item{\code{start}}{a single dissimilarity, or a list of
      dissimilarities to be employed as starting values.}
  }

  The default optimization using conjugate gradients should work
  reasonably well for medium to large size problems.  For \dQuote{small}
  ones, using \code{nlm} is usually faster.  Note that the number of
  ultrametric constraints is of the order \eqn{n^3}, where \eqn{n} is
  the number of objects in the dissimilarity object, suggesting to use
  the \acronym{SUMT} approach in favor of
  \code{\link[stats]{constrOptim}}.

  If starting values for the \acronym{SUMT} are provided via
  \code{start}, the number of starting values gives the number of runs
  to be performed, and control option \code{nruns} is ignored.
  Otherwise, \code{nruns} starting values are obtained by random shaking
  of the dissimilarity to be fitted.  In the case of multiple
  \acronym{SUMT} runs, the (first) best solution found is returned.

  Method \code{"IP"} implements the Iterative Projection approach of
  \bibcitet{Hubert+Arabie:1995}.  This iteratively projects the current
  dissimilarities to the closed convex set given by the ultrametric
  constraints (3-point conditions) for a single index triple \eqn{(i, j,
    k)}, in fact replacing the two largest values among \eqn{d_{ij},
    d_{ik}, d_{jk}} by their mean.  The following control parameters can
  be provided via the \code{control} argument.

  \describe{
    \item{\code{nruns}}{an integer giving the number of runs to be
      performed.  Defaults to 1.}
    \item{\code{order}}{a permutation of the numbers from 1 to the
      number of objects in \code{x}, specifying the order in which the
      ultrametric constraints are considered, or a list of such
      permutations.}
    \item{\code{maxiter}}{an integer giving the maximal number of
      iterations to be employed.}
    \item{\code{tol}}{a double indicating the maximal convergence
      tolerance.  The algorithm stops if the total absolute change in
      the dissimilarities in an iteration is less than \code{tol}.}
    \item{\code{verbose}}{a logical indicating whether to provide some
      output on minimization progress.  Defaults to
      \code{getOption("verbose")}.}
  }

  If permutations are provided via \code{order}, the number of these
  gives the number of runs to be performed, and control option
  \code{nruns} is ignored.  Otherwise, \code{nruns} randomly generated
  orders are tried.  In the case of multiple runs, the (first) best
  solution found is returned.

  Non-identical weights and incomplete dissimilarities are currently not
  supported.

  Method \code{"IR"} implements the Iterative Reduction approach
  suggested by \bibcitet{Roux:1988}, see also
  \bibcitet{Barthélemy+Guénoche:1991}.
  This is similar to the Iterative Projection method, but modifies the
  dissimilarities between objects proportionally to the aggregated
  change incurred from the ultrametric projections.  Available control
  parameters are identical to those of method \code{"IP"}.

  Non-identical weights and incomplete dissimilarities are currently not
  supported.

  It should be noted that all methods are heuristics which can not be
  guaranteed to find the global minimum.  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.
}
\references{
  \bibinfo{De_Soete:1984a}{year}{1984a}
  \bibinfo{De_Soete:1984b}{year}{1984b}  
  \bibshow{Barthélemy+Guénoche:1991,
    Carroll+Pruzansky:1980,
    Fiacco+McCormick:1968,
    Hubert+Arabie:1995,
    Křivánek+Morávek:1986,
    Roux:1988,
    De_Soete:1984a,
    De_Soete:1984b}
}
\seealso{
  \code{\link{cl_consensus}} for computing least squares (Euclidean)
  consensus hierarchies by least squares fitting of average ultrametric
  distances;
  \code{\link{l1_fit_ultrametric}}.
}
\examples{
## Least squares fit of an ultrametric to the Miller-Nicely consonant
## phoneme confusion data.
data("Phonemes")
## Note that the Phonemes data set has the consonant misclassification
## probabilities, i.e., the similarities between the phonemes.
d <- as.dist(1 - Phonemes)
u <- ls_fit_ultrametric(d, control = list(verbose = TRUE))
## Cophenetic correlation:
cor(d, u)
## Plot:
plot(u)
## ("Basically" the same as Figure 1 in de Soete (1984b).)
}
\keyword{cluster}
\keyword{optimize}
