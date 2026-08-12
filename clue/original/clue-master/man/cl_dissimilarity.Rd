\name{cl_dissimilarity}
\encoding{UTF-8}
\alias{cl_dissimilarity}
\title{Dissimilarity Between Partitions or Hierarchies}
\description{Compute the dissimilarity between (ensembles) of partitions
  or hierarchies.}
\usage{
cl_dissimilarity(x, y = NULL, method = "euclidean", \dots)
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies and dissimilarities,
    or something coercible to that (see \code{\link{cl_ensemble}}).}
  \item{y}{\code{NULL} (default), or as for \code{x}.}
  \item{method}{a character string specifying one of the built-in
    methods for computing dissimilarity, or a function to be taken as
    a user-defined method.  If a character string, its lower-cased
    version is matched against the lower-cased names of the available
    built-in methods using \code{\link{pmatch}}.  See \bold{Details} for
    available built-in methods.}
  \item{\dots}{further arguments to be passed to methods.}
}
\value{
  If \code{y} is \code{NULL}, an object of class
  \code{"cl_dissimilarity"} containing the dissimilarities between all
  pairs of components of \code{x}.  Otherwise, an object of class
  \code{"cl_cross_dissimilarity"} with the dissimilarities between the
  components of \code{x} and the components of \code{y}.
}
\details{
  If \code{y} is given, its components must be of the same kind as those
  of \code{x} (i.e., components must either all be partitions, or all be
  hierarchies or dissimilarities).

  If all components are partitions, the following built-in methods for
  measuring dissimilarity between two partitions with respective
  membership matrices \eqn{u} and \eqn{v} (brought to a common number of
  columns) are available:

  \describe{
    \item{\code{"euclidean"}}{the Euclidean dissimilarity of the
      memberships, i.e., the square root of the minimal sum of the
      squared differences of \eqn{u} and all column permutations of
      \eqn{v}.  See \bibcitet{Dimitriadou+Weingessel+Hornik:2002}.}
    \item{\code{"manhattan"}}{the Manhattan dissimilarity of the
      memberships, i.e., the minimal sum of the absolute differences of
      \eqn{u} and all column permutations of \eqn{v}.}
    \item{\code{"comemberships"}}{the Euclidean dissimilarity of the
      elements of the co-membership matrices \eqn{C(u) = u u'} and
      \eqn{C(v)}, i.e., the square root of the sum of the squared
      differences of \eqn{C(u)} and \eqn{C(v)}.}
    \item{\code{"symdiff"}}{the cardinality of the symmetric set
      difference of the sets of co-classified pairs of distinct objects
      in the partitions.  I.e., the number of distinct pairs of objects
      in the same class in exactly one of the partitions.
      (Alternatively, the cardinality of the symmetric set difference
      between the (binary) equivalence relations corresponding to the
      partitions.)  For soft partitions, (currently) the symmetric set
      difference of the corresponding nearest hard partitions is used.}
    \item{\code{"Rand"}}{the Rand distance, i.e., the rate of distinct
      pairs of objects in the same class in exactly one of the
      partitions.  (Related to the Rand index \eqn{a} via the linear
      transformation \eqn{d = (1 - a) / 2}.)  For soft partitions,
      (currently) the Rand distance of the corresponding nearest hard
      partitions is used.}
    \item{\code{"GV1"}}{the square root of the dissimilarity
      \eqn{\Delta_1}{Delta_1} used for the first model in
      \bibcitet{Gordon+Vichi:2001},
      i.e., the square root of the minimal sum of the
      squared differences of the \emph{matched} non-zero columns of
      \eqn{u} and \eqn{v}.}
    \item{\code{"BA/\var{d}"}}{distance measures for hard partitions
      discussed in \bibcitet{Boorman+Arabie:1972}, with \var{d} one of
      \samp{A}, \samp{C}, \samp{D}, or \samp{E}.  For soft partitions,
      the distances of the corresponding nearest hard partitions are
      used.

      \code{"BA/A"} is the minimum number of single element moves (move
      from one class to another or a new one) needed to transform one
      partition into the other.  Introduced in \bibcitet{Rubin:1967}.

      \code{"BA/C"} is the minimum number of lattice moves for
      transforming one partition into the other, where partitions are
      said to be connected by a lattice move if one is \emph{just} finer
      than the other (i.e., there is no other partition between them) in
      the partition lattice (see \code{\link{cl_meet}}).  Equivalently,
      with \eqn{z} the join of \code{x} and \code{y} and \eqn{S} giving
      the number of classes, this can be written as \eqn{S(x) + S(y) - 2
	S(z)}.  Attributed to David Pavy.
	
      \code{"BA/D"} is the \dQuote{pair-bonds} distance, which can be
      defined as \eqn{S(x) + S(y) - 2 S(z)}, with \eqn{z} the meet of
      \code{x} and \code{y} and \eqn{S} the \emph{supervaluation} (i.e.,
      non-decreasing with respect to the partial order on the partition
      lattice) function \eqn{\sum_i (n_i (n_i - 1)) / (n (n - 1))},
      where the \eqn{n_i} are the numbers of objects in the respective
      classes of the partition (such that \eqn{n_i (n_i - 1) / 2} are the
      numbers of pair bonds in the classes), and \eqn{n} the total
      number of objects.

      \code{"BA/E"} is the normalized information distance, defined as
      \eqn{1 - I / H}, where \eqn{I} is the average mutual information
      between the partitions, and \eqn{H} is the average entropy of the
      meet \eqn{z} of the partitions.
      Introduced in \bibcitet{Rajski:1961}.

      (Boorman and Arabie also discuss a distance measure (\eqn{B})
      based on the minimum number of set moves needed to transform one
      partition into the other, which, differently from the \eqn{A} and
      \eqn{C} distance measures is hard to compute \bibcitep{Day:1981}
      and (currently) not provided.)}
    \item{\code{"VI"}}{Variation of Information, see
      \bibcitet{Meila:2003}.  If
      \code{\dots} has an argument named \code{weights}, it is taken to
      specify case weights.}
    \item{\code{"Mallows"}}{the Mallows-type distance by
      \bibcitet{Zhou+Li+Zha:2005},
      which is related to the Monge-Kantorovich mass
      transfer problem, and given as the \eqn{p}-th root of the minimal
      value of the transportation problem \eqn{\sum w_{jk} \sum_i
      |u_{ij} - v_{ik}| ^ p} with constraints \eqn{w_{jk} \ge 0},
      \eqn{\sum_j w_{jk} = \alpha_j}, \eqn{\sum_k w_{jk} = \beta_k},
      where \eqn{\sum_j \alpha_j = \sum_k \beta_k}.  The parameters
      \eqn{p}, \eqn{\alpha} and \eqn{\beta} all default to one (in this
      case, the Mallows distance coincides with the Manhattan
      dissimilarity), and can be specified via additional arguments
      named \code{p}, \code{alpha}, and \code{beta}, respectively.}
    \item{\code{"CSSD"}}{the Cluster Similarity Sensitive Distance of
      \bibcitet{Zhou+Li+Zha:2005},
      which is given as the minimal value of
      \eqn{\sum_{k,l} (1 - 2 w_{kl} / (\alpha_k + \beta_l)) L_{kl}},
      where \eqn{L_{kl} = \sum_i u_{ik} v_{il} d(p_{x;k}, p_{y;l})} with
      \eqn{p_{x;k}} and \eqn{p_{y;l}} the prototype of the \eqn{k}-th
      class of \code{x} and the \eqn{l}-th class of \code{y},
      respectively, \eqn{d} is the distance between these, and the
      \eqn{w_{kl}} as for Mallows distance.  If prototypes are matrices,
      the Euclidean distance between these is used as default.  Using
      the additional argument \code{L}, one can give a matrix of
      \eqn{L_{kl}} values, or the function \eqn{d}.  Parameters
      \eqn{\alpha} and \eqn{\beta} all default to one, and can be
      specified via additional arguments named \code{alpha} and
      \code{beta}, respectively.}
  }

  For hard partitions, both Manhattan and squared Euclidean
  dissimilarity give twice the \emph{transfer distance}
  \bibcitep{Charon+Denoeud+Guénoche:2006},
  which is the minimum number of objects that must be removed so
  that the implied partitions (restrictions to the remaining objects)
  are identical.  This is also known as the \emph{\eqn{R}-metric} in
  \bibcitet{Day:1981},
  i.e., the number of augmentations and removals of single
  objects needed to transform one partition into the other, and the
  \emph{partition-distance} in \bibcitet{Gusfield:2002}, and equals
  twice the number of single element moves distance of Boorman and
  Arabie. 

  For hard partitions, the pair-bonds (Boorman-Arabie \eqn{D}) distance
  is identical to the Rand distance, and can also be written as the
  Manhattan distance between the co-membership matrices corresponding to
  the partitions, or equivalently, their symdiff distance, normalized by
  \eqn{n (n - 1)}.

  If all components are hierarchies, available built-in methods for
  measuring dissimilarity between two hierarchies with respective
  ultrametrics \eqn{u} and \eqn{v} are as follows.

  \describe{
    \item{\code{"euclidean"}}{the Euclidean dissimilarity of the
      ultrametrics (i.e., the square root of the sum of the squared
      differences of \eqn{u} and \eqn{v}).}
    \item{\code{"manhattan"}}{the Manhattan dissimilarity of the
      ultrametrics (i.e., the sum of the absolute differences of \eqn{u}
      and \eqn{v}).}
    \item{\code{"cophenetic"}}{\eqn{1 - c^2}, where \eqn{c} is the
      cophenetic correlation coefficient (i.e., the product-moment
      correlation of the ultrametrics).}
    \item{\code{"gamma"}}{the rate of inversions between the
      ultrametrics (i.e., the rate of pairs \eqn{(i,j)} and \eqn{(k,l)}
      for which \eqn{u_{ij} < u_{kl}} and \eqn{v_{ij} > v_{kl}}).}
    \item{\code{"symdiff"}}{the cardinality of the symmetric set
      difference of the sets of classes (hierarchies in the strict
      sense) induced by the dendrograms.  I.e., the number of sets of
      objects obtained by a split in exactly one of the hierarchies.}
    \item{\code{"Chebyshev"}}{the Chebyshev (maximal) dissimilarity of
      the ultrametrics (i.e., the maximum of the absolute differences of
      \eqn{u} and \eqn{v}).}
    \item{\code{"Lyapunov"}}{the logarithm of the product of the
      maximal and minimal ratios of the ultrametrics.  This is also
      known as the \dQuote{Hilbert projective metric} on the cone
      represented by the ultrametrics
      \bibcitep{e.g.,|Jardine+Sibson:1971|page 107},
      and only defined for \emph{strict} ultrametrics (which
      are strictly positive for distinct objects).}
    \item{\code{"BO"}}{the \eqn{m_\delta} family of tree metrics by
      \bibcitet{Boorman+Olivier:1973}, which are of the form
      \eqn{m_\delta = \int_0^\infty \delta(p(h), q(h)) dh}, where
      \eqn{p(h)} and
      \eqn{q(h)} are the hard partitions obtaining by cutting the trees
      (dendrograms) at height \eqn{h}, and \eqn{\delta} is a suitably
      dissimilarity measure for partitions.  In particular, when taking
      \eqn{\delta} as symdiff or Rand dissimilarity, \eqn{m_\delta} is
      the Manhattan dissimilarity of the hierarchies.

      If \code{\dots} has an argument named \code{delta} it is taken to
      specify the partition dissimilarity \eqn{\delta} to be employed.}
    \item{\code{"spectral"}}{the spectral norm (2-norm) of the
      differences of the ultrametrics, suggested in
      \bibcitet{Mérigot+Durbec+Gaertner:2010}.}
  }

  The measures based on ultrametrics also allow computing dissimilarity
  with \dQuote{raw} dissimilarities on the underlying objects (R objects
  inheriting from class \code{"dist"}).

  If a user-defined dissimilarity method is to be employed, it must be a
  function taking two clusterings as its arguments.

  Symmetric dissimilarity objects of class \code{"cl_dissimilarity"} are
  implemented as symmetric proximity objects with self-proximities
  identical to zero, and inherit from class \code{"cl_proximity"}.  They
  can be coerced to dense square matrices using \code{as.matrix}.  It
  is possible to use 2-index matrix-style subscripting for such objects;
  unless this uses identical row and column indices, this results in a
  (non-symmetric dissimilarity) object of class
  \code{"cl_cross_dissimilarity"}.

  Symmetric dissimilarity objects also inherit from class
  \code{"\link{dist}"} (although they currently do not \dQuote{strictly}
  extend this class), thus making it possible to use them directly for
  clustering algorithms based on dissimilarity matrices of this class,
  see the examples.
}
\references{
  \bibshow{Boorman+Arabie:1972,
    Boorman+Olivier:1973,
    Charon+Denoeud+Guénoche:2006,
    Day:1981,
    Dimitriadou+Weingessel+Hornik:2002,
    Gordon+Vichi:2001,
    Gusfield:2002,
    Jardine+Sibson:1971,
    Meila:2003,
    Mérigot+Durbec+Gaertner:2010,
    Rajski:1961,
    Rubin:1967,
    Zhou+Li+Zha:2005}
}  
\seealso{
  \code{\link{cl_agreement}}
}
\examples{
## An ensemble of partitions.
data("CKME")
pens <- CKME[1 : 30]
diss <- cl_dissimilarity(pens)
summary(c(diss))
cl_dissimilarity(pens[1:5], pens[6:7])
## Equivalently, using subscripting.
diss[1:5, 6:7]
## Can use the dissimilarities for "secondary" clustering
## (e.g. obtaining hierarchies of partitions):
hc <- hclust(diss)
plot(hc)

## Example from Boorman and Arabie (1972).
P1 <- as.cl_partition(c(1, 2, 2, 2, 3, 3, 2, 2))
P2 <- as.cl_partition(c(1, 1, 2, 2, 3, 3, 4, 4))
cl_dissimilarity(P1, P2, "BA/A")
cl_dissimilarity(P1, P2, "BA/C")

## Hierarchical clustering.
d <- dist(USArrests)
x <- hclust(d)
cl_dissimilarity(x, d, "cophenetic")
cl_dissimilarity(x, d, "gamma")
}
\keyword{cluster}
