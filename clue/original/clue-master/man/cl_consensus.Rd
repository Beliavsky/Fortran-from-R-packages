\name{cl_consensus}
\alias{cl_consensus}
\title{Consensus Partitions and Hierarchies}
\description{
  Compute the consensus clustering of an ensemble of partitions or
  hierarchies.
}
\usage{
cl_consensus(x, method = NULL, weights = 1, control = list())
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies, or something
    coercible to that (see \code{\link{cl_ensemble}}).}
  \item{method}{a character string specifying one of the built-in
    methods for computing consensus clusterings, or a function to be
    taken as a user-defined method, or \code{NULL} (default value).  If
    a character string, its lower-cased version is matched against the
    lower-cased names of the available built-in methods using
    \code{\link{pmatch}}.  See \bold{Details} for available built-in
    methods and defaults.}
  \item{weights}{a numeric vector with non-negative case weights.
    Recycled to the number of elements in the ensemble given by \code{x}
    if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  The consensus partition or hierarchy.
}
\details{
  Consensus clusterings \dQuote{synthesize} the information in the
  elements of a cluster ensemble into a single clustering, often by
  minimizing a criterion function measuring how dissimilar consensus
  candidates are from the (elements of) the ensemble (the so-called
  \dQuote{optimization approach} to consensus clustering).  

  The most popular criterion functions are of the form \eqn{L(x) = \sum
    w_b d(x_b, x)^p}, where \eqn{d} is a suitable dissimilarity measure
  (see \code{\link{cl_dissimilarity}}), \eqn{w_b} is the case weight
  given to element \eqn{x_b} of the ensemble, and \eqn{p \ge 1}.  If
  \eqn{p = 1} and minimization is over all possible base clusterings, a
  consensus solution is called a \emph{median} of the ensemble; if
  minimization is restricted to the elements of the ensemble, a
  consensus solution is called a \emph{medoid} (see
  \code{\link{cl_medoid}}).  For \eqn{p = 2}, we obtain \emph{least
    squares} consensus partitions and hierarchies (generalized means).
  See also \bibcitet{R:Gordon:1999} for more information.

  If all elements of the ensemble are partitions, the built-in consensus
  methods compute consensus partitions by minimizing a criterion of the
  form \eqn{L(x) = \sum w_b d(x_b, x)^p} over all hard or soft
  partitions \eqn{x} with a given (maximal) number \eqn{k} of classes.
  Available built-in methods are as follows.
  
  \describe{
    \item{\code{"SE"}}{a fixed-point algorithm for obtaining \emph{soft}
      least squares Euclidean consensus partitions (i.e., for minimizing
      \eqn{L} with Euclidean dissimilarity \eqn{d} and \eqn{p =	2} over
      all soft partitions with a given maximal number of classes).

      This iterates between individually matching all partitions to the
      current approximation to the consensus partition, and computing
      the next approximation as the membership matrix closest to a
      suitable weighted average of the memberships of all partitions
      after permuting their columns for the optimal matchings of class
      ids.

      The following control parameters are available for this method.
      
      \describe{
	\item{\code{k}}{an integer giving the number of classes to be
	  used for the least squares consensus partition.
	  By default, the maximal number of classes in the ensemble is
	  used.}
	\item{\code{maxiter}}{an integer giving the maximal number of
	  iterations to be performed.
	  Defaults to 100.}
	\item{\code{nruns}}{an integer giving the number of runs to be
	  performed.  Defaults to 1.}
	\item{\code{reltol}}{the relative convergence tolerance.
	  Defaults to \code{sqrt(.Machine$double.eps)}.}
	\item{\code{start}}{a matrix with number of rows equal to the
	  number of objects of the cluster ensemble, and \eqn{k}
	  columns, to be used as a starting value, or a list of such
	  matrices.  By default, suitable random membership matrices are
	  used.}
	\item{\code{verbose}}{a logical indicating whether to provide
	  some output on minimization progress.
	  Defaults to \code{getOption("verbose")}.}
      }
      In the case of multiple runs, the first optimum found is returned.

      This method can also be referred to as \code{"soft/euclidean"}.
    }

    \item{\code{"GV1"}}{the fixed-point algorithm for the \dQuote{first
	model} in \bibcitet{Gordon+Vichi:2001} for minimizing \eqn{L}
      with \eqn{d} being GV1 dissimilarity and \eqn{p = 2} over all soft
      partitions with a given maximal number of classes.

      This is similar to \code{"SE"}, but uses GV1 rather than Euclidean
      dissimilarity.

      Available control parameters are the same as for \code{"SE"}.
    }
    
    \item{\code{"DWH"}}{an extension of the greedy algorithm in
      \bibcitet{Dimitriadou+Weingessel+Hornik:2002} for (approximately)
      obtaining soft least squares Euclidean consensus partitions.
      The reference provides some structure theory relating finding
      the consensus partition to an instance of the multiple assignment
      problem, which is known to be NP-hard, and suggests a simple
      heuristic based on successively matching an individual partition
      \eqn{x_b} to the current approximation to the consensus partition,
      and compute the memberships of the next approximation as a
      weighted average of those of the current one and of \eqn{x_b}
      after permuting its columns for the optimal matching of class
      ids.

      The following control parameters are available for this method.
      
      \describe{
	\item{\code{k}}{an integer giving the number of classes to be
	  used for the least squares consensus partition.  By default,
	  the maximal number of classes in the ensemble is used.}
	\item{\code{order}}{a permutation of the integers from 1 to the
	  size of the ensemble, specifying the order in which the
	  partitions in the ensemble should be aggregated.  Defaults to
	  using a random permutation (unlike the reference, which does
	  not permute at all).}
      }
    }

    \item{\code{"HE"}}{a fixed-point algorithm for obtaining \emph{hard}
      least squares Euclidean consensus partitions (i.e., for minimizing
      \eqn{L} with Euclidean dissimilarity \eqn{d} and \eqn{p =	2} over
      all hard partitions with a given maximal number of classes.)

      Available control parameters are the same as for \code{"SE"}.

      This method can also be referred to as \code{"hard/euclidean"}.
    }
    
    \item{\code{"SM"}}{a fixed-point algorithm for obtaining \emph{soft}
      median Manhattan consensus partitions (i.e., for minimizing
      \eqn{L} with Manhattan dissimilarity \eqn{d} and \eqn{p =	1} over 
      all soft partitions with a given maximal number of classes).

      Available control parameters are the same as for \code{"SE"}.

      This method can also be referred to as \code{"soft/manhattan"}.
    }

    \item{\code{"HM"}}{a fixed-point algorithm for obtaining \emph{hard}
      median Manhattan consensus partitions (i.e., for minimizing
      \eqn{L} with Manhattan dissimilarity \eqn{d} and \eqn{p =	1} over 
      all hard partitions with a given maximal number of classes).

      Available control parameters are the same as for \code{"SE"}.

      This method can also be referred to as \code{"hard/manhattan"}.
    }

    \item{\code{"GV3"}}{a \acronym{SUMT} algorithm for the \dQuote{third
      model} in \bibcitet{Gordon+Vichi:2001} for minimizing \eqn{L} with
      \eqn{d} being co-membership dissimilarity and \eqn{p = 2}.  (See
      \code{\link{sumt}} for more information on the \acronym{SUMT}
      approach.)  This optimization problem is equivalent to finding the
      membership matrix \eqn{m} for which the sum of the squared
      differences between \eqn{C(m) = m m'} and the weighted average
      co-membership matrix \eqn{\sum_b w_b C(m_b)} of the partitions is
      minimal.

      Available control parameters are \code{method}, \code{control},
      \code{eps}, \code{q}, and \code{verbose}, which have the same
      roles as for \code{\link{sumt}}, and the following.
      
      \describe{
	\item{\code{k}}{an integer giving the number of classes to be
	  used for the least squares consensus partition.  By default,
	  the maximal number of classes in the ensemble is used.}
	\item{\code{nruns}}{an integer giving the number of runs to be
	  performed.  Defaults to 1.}
	\item{\code{start}}{a matrix with number of rows equal to the
	  size of the cluster ensemble, and \eqn{k} columns, to be used
	  as a starting value, or a list of such matrices.  By default,
	  a membership based on a rank \eqn{k} approximation to the
	  weighted average co-membership matrix is used.}
      }
      In the case of multiple runs, the first optimum found is returned.
      
    }

    \item{\code{"soft/symdiff"}}{a \acronym{SUMT} approach for
      minimizing \eqn{L = \sum w_b d(x_b, x)} over all soft partitions
      with a given maximal number of classes, where \eqn{d} is the
      Manhattan dissimilarity of the co-membership matrices (coinciding
      with symdiff partition dissimilarity in the case of hard
      partitions).

      Available control parameters are the same as for \code{"GV3"}.
    }

    \item{\code{"hard/symdiff"}}{an exact solver for minimizing
      \eqn{L = \sum w_b d(x_b, x)} over all hard partitions (possibly
      with a given maximal number of classes as specified by the control
      parameter \code{k}), where \eqn{d} is symdiff partition
      dissimilarity (so that soft partitions in the ensemble are
      replaced by their closest hard partitions), or equivalently, Rand
      distance or pair-bonds (Boorman-Arabie \eqn{D}) distance.  The
      consensus solution is found via mixed linear or quadratic
      programming.
    }
  }
  
  By default, method \code{"SE"} is used for ensembles of partitions.
  
  If all elements of the ensemble are hierarchies, the following
  built-in methods for computing consensus hierarchies are available.
  
  \describe{
    \item{\code{"euclidean"}}{an algorithm for minimizing
      \eqn{L(x) = \sum w_b d(x_b, x) ^ 2} over all dendrograms, where
      \eqn{d} is Euclidean dissimilarity.  This is equivalent to finding
      the best least squares ultrametric approximation of the weighted
      average \eqn{d = \sum w_b u_b} of the ultrametrics \eqn{u_b} of
      the hierarchies \eqn{x_b}, which is attempted by calling
      \code{\link{ls_fit_ultrametric}} on \eqn{d} with appropriate
      control parameters.

      This method can also be referred to as \code{"cophenetic"}.
    }

    \item{\code{"manhattan"}}{a \acronym{SUMT} for minimizing
      \eqn{L = \sum w_b d(x_b, x)} over all dendrograms, where \eqn{d}
      is Manhattan dissimilarity.

      Available control parameters are the same as for
      \code{"euclidean"}.
    }
    
    \item{\code{"majority"}}{a hierarchy obtained from an extension of
      the majority consensus tree of \bibcitet{Margush+McMorris:1981},
      which minimizes \eqn{L(x) = \sum w_b d(x_b, x)} over all dendrograms,
      where \eqn{d} is the symmetric difference dissimilarity.  The
      unweighted \eqn{p}-majority tree is the \eqn{n}-tree (hierarchy in
      the strict sense) consisting of all subsets of objects contained
      in more than \eqn{100 p} percent of the \eqn{n}-trees \eqn{T_b}
      induced by the dendrograms, where \eqn{1/2 \le p < 1} and
      \eqn{p = 1/2} (default) corresponds to the standard majority tree.
      In the weighted case, it consists of all subsets \eqn{A} for which
      \eqn{\sum_{b: A \in T_b} w_b > W p}, where \eqn{W = \sum_b w_b}.
      We also allow for \eqn{p = 1}, which gives the \emph{strict
	consensus tree} consisting of all subsets contained in each of
      the \eqn{n}-trees.  The majority dendrogram returned is a
      representation of the majority tree where all splits have height
      one.

      The fraction \eqn{p} can be specified via the control parameter
      \code{p}.
    }
  }

  By default, method \code{"euclidean"} is used for ensembles of
  hierarchies.

  If a user-defined consensus method is to be employed, it must be a
  function taking the cluster ensemble, the case weights, and a list of
  control parameters as its arguments, with formals named \code{x},
  \code{weights}, and \code{control}, respectively.

  Most built-in methods use heuristics for solving hard optimization
  problems, and cannot be guaranteed to find a global minimum.  Standard
  practice would recommend to use the best solution found in
  \dQuote{sufficiently many} replications of the methods.
}
\references{
  \bibshow{Dimitriadou+Weingessel+Hornik:2002,
    Gordon+Vichi:2001,
    R:Gordon:1999,
    Margush+McMorris:1981}
}
\seealso{
  \code{\link{cl_medoid}},
  \code{\link[ape]{consensus}}
}
\examples{
## Consensus partition for the Rosenberg-Kim kinship terms partition
## data based on co-membership dissimilarities.
data("Kinship82")
m1 <- cl_consensus(Kinship82, method = "GV3",
                   control = list(k = 3, verbose = TRUE))
## (Note that one should really use several replicates of this.)
## Value for criterion function to be minimized:
sum(cl_dissimilarity(Kinship82, m1, "comem") ^ 2)
## Compare to the consensus solution given in Gordon & Vichi (2001).
data("Kinship82_Consensus")
m2 <- Kinship82_Consensus[["JMF"]]
sum(cl_dissimilarity(Kinship82, m2, "comem") ^ 2)
## Seems we get a better solution ...
## How dissimilar are these solutions?
cl_dissimilarity(m1, m2, "comem")
## How "fuzzy" are they?
cl_fuzziness(cl_ensemble(m1, m2))
## Do the "nearest" hard partitions fully agree?
cl_dissimilarity(as.cl_hard_partition(m1),
                 as.cl_hard_partition(m2))

## Consensus partition for the Gordon and Vichi (2001) macroeconomic
## partition data based on Euclidean dissimilarities.
data("GVME")
set.seed(1)
## First, using k = 2 classes.
m1 <- cl_consensus(GVME, method = "GV1",
                   control = list(k = 2, verbose = TRUE))
## (Note that one should really use several replicates of this.)
## Value of criterion function to be minimized:
sum(cl_dissimilarity(GVME, m1, "GV1") ^ 2)
## Compare to the consensus solution given in Gordon & Vichi (2001).
data("GVME_Consensus")
m2 <- GVME_Consensus[["MF1/2"]]
sum(cl_dissimilarity(GVME, m2, "GV1") ^ 2)
## Seems we get a slightly  better solution ...
## But note that
cl_dissimilarity(m1, m2, "GV1")
## and that the maximal deviation of the memberships is
max(abs(cl_membership(m1) - cl_membership(m2)))
## so the differences seem to be due to rounding.
## Do the "nearest" hard partitions fully agree?
table(cl_class_ids(m1), cl_class_ids(m2))

## And now for k = 3 classes.
m1 <- cl_consensus(GVME, method = "GV1",
                   control = list(k = 3, verbose = TRUE))
sum(cl_dissimilarity(GVME, m1, "GV1") ^ 2)
## Compare to the consensus solution given in Gordon & Vichi (2001).
m2 <- GVME_Consensus[["MF1/3"]]
sum(cl_dissimilarity(GVME, m2, "GV1") ^ 2)
## This time we look much better ...
## How dissimilar are these solutions?
cl_dissimilarity(m1, m2, "GV1")
## Do the "nearest" hard partitions fully agree?
table(cl_class_ids(m1), cl_class_ids(m2))
}
\keyword{cluster}
