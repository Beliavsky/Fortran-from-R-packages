\name{pclust}
\alias{pclust}
\alias{pclust_family}
\alias{pclust_object}
\title{Prototype-Based Partitioning}
\description{
  Obtain prototype-based partitions of objects by minimizing the criterion
  \eqn{\sum w_b u_{bj}^m d(x_b, p_j)^e}, the sum of the case-weighted and
  membership-weighted \eqn{e}-th powers of the dissimilarities between
  the objects \eqn{x_b} and the prototypes \eqn{p_j}, for suitable
  dissimilarities \eqn{d} and exponents \eqn{e}.
}
\usage{
pclust(x, k, family, m = 1, weights = 1, control = list())
pclust_family(D, C, init = NULL, description = NULL, e = 1,
              .modify = NULL, .subset = NULL)
pclust_object(prototypes, membership, cluster, family, m = 1,
              value, ..., classes = NULL, attributes = NULL)
}
\arguments{
  \item{x}{the object to be partitioned.}
  \item{k}{an integer giving the number of classes to be used in the
    partition.}
  \item{family}{an object of class \code{"pclust_family"} as generated
    by \code{pclust_family}, containing the information about \eqn{d}
    and \eqn{e}.}
  \item{m}{a number not less than 1 controlling the softness of the
    partition (as the \dQuote{fuzzification parameter} of the fuzzy
    \eqn{c}-means algorithm).  The default value of 1 corresponds to
    hard partitions obtained from a generalized \eqn{k}-means problem;
    values greater than one give partitions of increasing softness
    obtained from a generalized fuzzy \eqn{c}-means problem.}
  \item{weights}{a numeric vector of non-negative case weights.
    Recycled to the number of elements given by \code{x} if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
  \item{D}{a function for computing dissimilarities \eqn{d} between
    objects and prototypes.}
  \item{C}{a \sQuote{consensus} function with formals \code{x},
    \code{weights} and \code{control} for computing a consensus
    prototype \eqn{p} minimizing \eqn{\sum_b w_b d(x_b, p) ^ e}.}
  \item{init}{a function with formals \code{x} and \code{k} initializing
    an object with \eqn{k} prototypes from the object \code{x} to be
    partitioned.}
  \item{description}{a character string describing the family.}
  \item{e}{a number giving the exponent \eqn{e} of the criterion.}
  \item{.modify}{a function with formals \code{x}, \code{i} and
    \code{value} for modifying a single prototype,
    or \code{NULL} (default).}
  \item{.subset}{a function with formals \code{x} and \code{i} for
    subsetting prototypes,
    or \code{NULL} (default).}
  \item{prototypes}{an object representing the prototypes of the
    partition.}
  \item{membership}{an object of class \code{"\link{cl_membership}"}
    with the membership values \eqn{u_{bj}}.}
  \item{cluster}{the class ids of the nearest hard partition.}
  \item{value}{the value of the criterion to be minimized.}
  \item{...}{further elements to be included in the generated pclust
    object.}
  \item{classes}{a character vector giving further classes to be given
    to the generated pclust object in addition to \code{"pclust"}, or
    \code{NULL} (default).}
  \item{attributes}{a list of attributes, or \code{NULL} (default).}
}
\value{
  \code{pclust} returns the partition found as an object of class
  \code{"pclust"} (as obtained by calling \code{pclust_object}) which in
  addition to the \emph{default} components contains \code{call} (the 
  matched call) and a \code{converged} attribute indicating convergence
  status (i.e., whether the maximal number of iterations was reached).

  \code{pclust_family} returns an object of class
  \code{"pclust_family"}, which is a list with components corresponding
  to the formals of \code{pclust_family}.

  \code{pclust_object} returns an object inheriting from class
  \code{"pclust"}, which is a list with components corresponding
  to the formals (up to and including \code{...}) of
  \code{pclust_object}, and additional classes and attributes specified
  by \code{classes} and \code{attributes}, respectively.
}
\details{
  For \eqn{m = 1}, a generalization of the Lloyd-Forgy variant of the
  \eqn{k}-means algorithm is used, which iterates between reclassifying
  objects to their closest prototypes (according to the dissimilarities
  given by \code{D}), and computing new prototypes as the consensus for
  the classes (using \code{C}).

  For \eqn{m > 1}, a generalization of the fuzzy \eqn{c}-means recipe
  \bibcitep{e.g.,|Bezdek:1981|} is used, which alternates between computing
  optimal memberships for fixed prototypes, and computing new prototypes
  as the suitably weighted consensus clusterings for the classes.

  This procedure is repeated until convergence occurs, or the maximal
  number of iterations is reached.

  Currently, no local improvement heuristics are provided.

  It is possible to perform several runs of the procedure via control
  arguments \code{nruns} or \code{start} (the default is to perform a
  single run), in which case the first partition with the smallest
  value of the criterion is returned.

  The dissimilarity and consensus functions as well as the exponent
  \eqn{e} are specified via \code{family}.  In principle, arbitrary
  representations of objects to be partitioned and prototypes (which do
  not necessarily have to be \dQuote{of the same kind}) can be employed.
  In addition to \code{D} and \code{C}, what is needed are means to
  obtain an initial collection of \eqn{k} prototypes (\code{init}), to
  modify a single prototype (\code{.modify}), and subset the prototypes
  (\code{.subset}).  By default, list and (currently, only dense) matrix
  (with the usual convention that the rows correspond to the objects)
  are supported.  Otherwise, the family has to provide the functions
  needed. 

  Available control parameters are as follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      iterations to be performed.
      Defaults to 100.}
    \item{\code{nruns}}{an integer giving the number of runs to be
      performed.
      Defaults to 1.} 
    \item{\code{reltol}}{the relative convergence tolerance.
      Defaults to \code{sqrt(.Machine$double.eps)}.}
    \item{\code{start}}{a list of prototype objects to be used as
      starting values.}
    \item{\code{verbose}}{a logical indicating whether to provide
      some output on minimization progress.
      Defaults to \code{getOption("verbose")}.}
    \item{\code{control}}{control parameters to be used in the consensus
      function.}
  }

  The fixed point approach employed is a heuristic which cannot be
  guaranteed to find the global minimum, in particular if \code{C} is
  not an exact minimizer.  Standard practice would recommend to use the
  best solution found in \dQuote{sufficiently many} replications of the
  base algorithm.
}
\references{
  \bibshow{Bezdek:1981}
}
\seealso{
  \code{\link[stats]{kmeans}},
  \code{\link[e1071]{cmeans}}.
}  
