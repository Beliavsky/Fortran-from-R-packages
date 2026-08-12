\name{cl_agreement}
\alias{cl_agreement}
\title{Agreement Between Partitions or Hierarchies}
\description{Compute the agreement between (ensembles) of partitions or
  hierarchies.
}
\usage{
cl_agreement(x, y = NULL, method = "euclidean", \dots)
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies and dissimilarities,
    or something coercible to that (see \code{\link{cl_ensemble}}).}
  \item{y}{\code{NULL} (default), or as for \code{x}.}
  \item{method}{a character string specifying one of the built-in
    methods for computing agreement, or a function to be taken as
    a user-defined method.  If a character string, its lower-cased
    version is matched against the lower-cased names of the available
    built-in methods using \code{\link{pmatch}}.  See \bold{Details} for
    available built-in methods.}
  \item{\dots}{further arguments to be passed to methods.}
}
\value{
  If \code{y} is \code{NULL}, an object of class \code{"cl_agreement"}
  containing the agreements between the all pairs of components of
  \code{x}.  Otherwise, an object of class \code{"cl_cross_agreement"}
  with the agreements between the components of \code{x} and the
  components of \code{y}.
}
\details{
  If \code{y} is given, its components must be of the same kind as those
  of \code{x} (i.e., components must either all be partitions, or all be
  hierarchies or dissimilarities).

  If all components are partitions, the following built-in methods for
  measuring agreement between two partitions with respective membership
  matrices \eqn{u} and \eqn{v} (brought to a common number of columns)
  are available:
  
  \describe{
    \item{\code{"euclidean"}}{\eqn{1 - d / m}, where \eqn{d} is the
      Euclidean dissimilarity of the memberships, i.e., the square root
      of the minimal sum of the squared differences of \eqn{u} and all
      column permutations of \eqn{v}, and \eqn{m} is an upper bound for
      the maximal Euclidean dissimilarity.
      See \bibcitet{Dimitriadou+Weingessel+Hornik:2002}.}
    \item{\code{"manhattan"}}{\eqn{1 - d / m}, where \eqn{d} is the
      Manhattan dissimilarity of the memberships, i.e., the minimal
      sum of the absolute differences of \eqn{u} and all column
      permutations of \eqn{v}, and \eqn{m} is an upper bound for the
      maximal Manhattan dissimilarity.}
    \item{\code{"Rand"}}{the Rand index (the rate of distinct pairs of
      objects both in the same class or both in different classes in
      both partitions), see \bibcitet{Rand:1971} or
      \bibcitet{|R:Gordon:1999|page 198}.
      For soft partitions, (currently) the Rand index of the
      corresponding nearest hard partitions is used.}
    \item{\code{"cRand"}}{the Rand index corrected for agreement by
      chance, see \bibcitet{Hubert+Arabie:1985} or
      \bibcitet{|R:Gordon:1999|page 198}.
      Can only be used for hard partitions.}
    \item{\code{"NMI"}}{Normalized Mutual Information, see
      \bibcitet{Strehl+Ghosh:2002}.
      For soft partitions, (currently) the NMI of the
      corresponding nearest hard partitions is used.} 
    \item{\code{"KP"}}{the Katz-Powell index, i.e., the product-moment
      correlation coefficient between the elements of the co-membership
      matrices \eqn{C(u) = u u'} and \eqn{C(v)}, respectively,
      see \bibcitet{Katz+Powell:1953}.  For soft partitions, (currently) the
      Katz-Powell index of the corresponding nearest hard partitions is
      used.  (Note that for hard partitions, the \eqn{(i,j)} entry of
      \eqn{C(u)} is one iff objects \eqn{i} and \eqn{j} are in the same
      class.)}
    \item{\code{"angle"}}{the maximal cosine of the angle between the
      elements of \eqn{u} and all column permutations of \eqn{v}.}
    \item{\code{"diag"}}{the maximal co-classification rate, i.e., the
      maximal rate of objects with the same class ids in both
      partitions after arbitrarily permuting the ids.}
    \item{\code{"FM"}}{the index of \bibcitet{Fowlkes+Mallows:1983}, i.e.,
      the ratio \eqn{N_{xy} / \sqrt{N_x N_y}}{N_xy / sqrt(N_x N_y)} of
      the number \eqn{N_{xy}}{N_xy} of distinct pairs of objects in the
      same class in both partitions and the geometric mean of the
      numbers \eqn{N_x} and \eqn{N_y} of distinct pairs of objects in
      the same class in partition \eqn{x} and partition \eqn{y},
      respectively.
      For soft partitions, (currently) the Fowlkes-Mallows index of the
      corresponding nearest hard partitions is used.}
    \item{\code{"Jaccard"}}{the Jaccard index, i.e., the ratio of the
      numbers of distinct pairs of objects in the same class in both
      partitions and in at least one partition, respectively.
      For soft partitions, (currently) the Jaccard index of the
      corresponding nearest hard partitions is used.}
    \item{\code{"purity"}}{the purity of the classes of \code{x} with
      respect to those of \code{y}, i.e.,
      \eqn{\sum_j \max_i n_{ij} / n},
      where \eqn{n_{ij}} is the joint frequency of objects in class
      \eqn{i} for \code{x} and in class \eqn{j} for \code{y}, and
      \eqn{n} is the total number of objects.}
    \item{\code{"PS"}}{Prediction Strength, see
      \bibcitet{Tibshirani+Walther:2005}: the minimum, over all classes
      \eqn{j} of \code{y}, of the maximal rate of objects in the same
      class for \code{x} and in class \eqn{j} for \code{y}.}
  }

  If all components are hierarchies, available built-in methods for
  measuring agreement between two hierarchies with respective
  ultrametrics \eqn{u} and \eqn{v} are as follows.

  \describe{
    \item{\code{"euclidean"}}{\eqn{1 / (1 + d)}, where \eqn{d} is the
      Euclidean dissimilarity of the ultrametrics (i.e., the square root
      of the sum of the squared differences of \eqn{u} and \eqn{v}).}
    \item{\code{"manhattan"}}{\eqn{1 / (1 + d)}, where \eqn{d} is the
      Manhattan dissimilarity of the ultrametrics (i.e., the sum of the
      absolute differences of \eqn{u} and \eqn{v}).}
    \item{\code{"cophenetic"}}{The cophenetic correlation coefficient.
      (I.e., the product-moment correlation of the ultrametrics.)}
    \item{\code{"angle"}}{the cosine of the angle between the
      ultrametrics.}
    \item{\code{"gamma"}}{\eqn{1 - d}, where \eqn{d} is the rate of
      inversions between the associated ultrametrics (i.e., the rate of
      pairs \eqn{(i,j)} and \eqn{(k,l)} for which \eqn{u_{ij} < u_{kl}}
      and \eqn{v_{ij} > v_{kl}}).  (This agreement measure is a linear
      transformation of Kruskal's \eqn{\gamma}{gamma}.)}
  }

  The measures based on ultrametrics also allow computing agreement with
  \dQuote{raw} dissimilarities on the underlying objects (R objects
  inheriting from class \code{"dist"}).

  If a user-defined agreement method is to be employed, it must be a
  function taking two clusterings as its arguments.

  Symmetric agreement objects of class \code{"cl_agreement"} are
  implemented as symmetric proximity objects with self-proximities
  identical to one, and inherit from class \code{"cl_proximity"}.  They
  can be coerced to dense square matrices using \code{as.matrix}.  It is
  possible to use 2-index matrix-style subscripting for such objects;
  unless this uses identical row and column indices, this results in a
  (non-symmetric agreement) object of class \code{"cl_cross_agreement"}.
}
\references{
  \bibshow{Dimitriadou+Weingessel+Hornik:2002,
    Fowlkes+Mallows:1983,
    R:Gordon:1999,
    Hubert+Arabie:1985,
    Katz+Powell:1953,
    Rand:1971,
    Strehl+Ghosh:2002,
    Tibshirani+Walther:2005}
}
\seealso{
  \code{\link{cl_dissimilarity}};
  \code{\link[e1071]{classAgreement}} in package \pkg{e1071}.
}
\examples{
## An ensemble of partitions.
data("CKME")
pens <- CKME[1 : 20]		# for saving precious time ...
summary(c(cl_agreement(pens)))
summary(c(cl_agreement(pens, method = "Rand")))
summary(c(cl_agreement(pens, method = "diag")))
cl_agreement(pens[1:5], pens[6:7], method = "NMI")
## Equivalently, using subscripting.
cl_agreement(pens, method = "NMI")[1:5, 6:7]

## An ensemble of hierarchies.
d <- dist(USArrests)
hclust_methods <-
    c("ward", "single", "complete", "average", "mcquitty")
hclust_results <- lapply(hclust_methods, function(m) hclust(d, m))
names(hclust_results) <- hclust_methods 
hens <- cl_ensemble(list = hclust_results)
summary(c(cl_agreement(hens)))
## Note that the Euclidean agreements are *very* small.
## This is because the ultrametrics differ substantially in height:
u <- lapply(hens, cl_ultrametric)
round(sapply(u, max), 3)
## Rescaling the ultrametrics to [0, 1] gives:
u <- lapply(u, function(x) (x - min(x)) / (max(x) - min(x)))
shens <- cl_ensemble(list = lapply(u, as.cl_dendrogram))
summary(c(cl_agreement(shens)))
## Au contraire ...
summary(c(cl_agreement(hens, method = "cophenetic")))
cl_agreement(hens[1:3], hens[4:5], method = "gamma")
}
\keyword{cluster}
