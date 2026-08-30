\name{dist}
\alias{dist}
\alias{simil}
\alias{print.simil}
\alias{print.dist}
\alias{print.crosssimil}
\alias{print.crossdist}
\alias{as.matrix.dist}
\alias{as.matrix.simil}
\alias{pr_simil2dist}
\alias{pr_dist2simil}
\alias{as.matrix}
\alias{as.dist}
\alias{as.simil}
\alias{row.dist}
\alias{col.dist}
\title{Matrix Distance/Similarity Computation}
\description{
  These functions compute and return the auto-distance/similarity matrix
  between either rows or columns of a matrix/data frame, or a list,
  as well as the cross-distance matrix between two matrices/data frames/lists.
}
\usage{
dist(x, y = NULL, method = NULL, ..., diag = FALSE, upper = FALSE,
     pairwise = FALSE, by_rows = TRUE, convert_similarities = TRUE,
     auto_convert_data_frames = TRUE)
simil(x, y = NULL, method = NULL, ..., diag = FALSE, upper = FALSE,
      pairwise = FALSE, by_rows = TRUE, convert_distances = TRUE,
      auto_convert_data_frames = TRUE)

pr_dist2simil(x)
pr_simil2dist(x)

as.dist(x, FUN = NULL)
as.simil(x, FUN = NULL)

\method{as.matrix}{dist}(x, diag = 0, \dots)
\method{as.matrix}{simil}(x, diag = NA, \dots)
}
\arguments{
  \item{x}{For \code{dist} and \code{simil}, a numeric matrix object, a data frame, or a list. A vector
    will be converted into a column matrix. For \code{as.simil} and
    \code{as.dist}, an object of class \code{dist} and
    \code{simil}, respectively, or a numeric matrix. For
    \code{pr_dist2simil} and \code{pr_simil2dist}, any numeric vector.}
  \item{y}{\code{NULL}, or a similar object than \code{x}}
  \item{method}{a function, a registry entry, or a mnemonic string referencing the
    proximity measure. A list of all available measures can be obtained
    using \code{\link{pr_DB}} (see examples). The default for \code{dist} is
    \code{"Euclidean"}, and for \code{simil} \code{"correlation"}.}
  \item{diag}{logical value indicating whether the diagonal of the
    distance/similarity matrix should be printed by
    \code{\link{print.dist}}/\code{\link{print.simil}}. Note that the
    diagonal values are never stored in \code{dist} objects. 

    In the context of \code{as.matrix} the value to use on the diagonal
    representing self-proximities. In case of similarities, this
    defaults to \code{NA} since a priori there are no upper bounds, so
    the maximum similarity needs to be specified by the user.}
  \item{upper}{logical value indicating whether the upper triangle of the
    distance/similarity matrix should be printed by
    \code{\link{print.dist}}/\code{\link{print.simil}}}
  \item{pairwise}{logical value indicating whether distances should be
    computed for the pairs of \code{x} and \code{y} only.}
  \item{by_rows}{logical indicating whether proximities between rows, or
    columns should be computed.}
  \item{convert_similarities, convert_distances}{logical indicating
    whether distances should be automatically converted into
    similarities (and the other way round) if needed.}
  \item{auto_convert_data_frames}{logical indicating whether data frames
    should be converted to matrices if all variables are numeric,
    or all are logical, or all are complex.}
  \item{FUN}{optional function to be used by \code{as.dist} and
    \code{as.simil}. If \code{NULL}, it is looked up in the method
    registry. If there is none specified there, \code{FUN} defaults to
    \code{pr_simil2dist} and \code{pr_dist2simil}, respectively.}
  \item{\dots}{further arguments passed to the proximity function.}
}
\details{
  The interface is fashioned after \code{\link[stats]{dist}}, but can
  also compute cross-distances, and allows user extensions by means of
  registry of all proximity measures (see \code{\link{pr_DB}}).

  Missing values are allowed but are excluded from all computations 
  involving the rows within which they occur. If some columns are
  excluded in calculating a Euclidean, Manhattan, Canberra or
  Minkowski distance, the sum is scaled up proportionally to the
  number of columns used (compare \code{\link[stats]{dist}} in
  package \pkg{stats}).

  Data frames are silently coerced to matrix if all columns are of
  (same) mode \code{numeric} or \code{logical}.

  Distance measures can be used with \code{simil}, and similarity
  measures with \code{dist}. In these cases, the result is transformed
  accordingly using the specified coercion functions (default:
  \eqn{\mathrm{pr\_simil2dist}(x) = 1 - \mathrm{abs}(x)}{pr_simil2dist(x) = 1 - abs(x)} and \eqn{\mathrm{pr\_dist2simil}(x) = 1 / (1 + x)}{pr_dist2simil(x) = 1 / (1 + x)}).
  Objects of class \code{simil} and \code{dist} can be converted one in
  another using \code{as.dist} and \code{as.simil}, respectively.

  Distance and similarity objects can conveniently be subset
  (see examples). Note that duplicate indexes are silently ignored.
}
\value{
  Auto distances/similarities are returned as an object of class \code{dist}/\code{simil} and 
  cross-distances/similarities as an object of class \code{crossdist}/\code{crosssimil}. 
}
\references{
  Anderberg, M.R. (1973), \emph{Cluster analysis for applications},
  359 pp., Academic Press, New York, NY, USA.
  
  Cox, M.F. and Cox, M.A.A. (2001), \emph{Multidimensional Scaling},
  Chapman and Hall.
  
  Sokol, R.S. and Sneath P.H.A (1963), \emph{Principles of Numerical
  Taxonomy}, W. H. Freeman and Co., San Francisco.
}
\author{David Meyer \email{David.Meyer@R-project.org}
  and Christian Buchta \email{Christian.Buchta@wu-wien.ac.at}}

\seealso{\code{\link[stats]{dist}} for compatibility information, and
  \code{\link{pr_DB}} for the proximity data base.}
\examples{
### show available proximities
summary(pr_DB)

### get more information about a particular one
pr_DB$get_entry("Jaccard")

### binary data
x <- matrix(sample(c(FALSE, TRUE), 8, replace = TRUE), ncol = 2)
dist(x, method = "Jaccard")

### for real-valued data
dist(x, method = "eJaccard")

### for positive real-valued data
dist(x, method = "fJaccard")

### cross distances
dist(x, x, method = "Jaccard")

### pairwise (diagonal)
dist(x, x, method = "Jaccard", 
	 pairwise = TRUE)

### this is the same but less efficient
as.matrix(stats::dist(x, method = "binary"))

### numeric data
x <- matrix(rnorm(16), ncol = 4)

## test inheritance of names
rownames(x) <- LETTERS[1:4]
colnames(x) <- letters[1:4]
dist(x)
dist(x, x)

## custom distance function
f <- function(x, y) sum(x * y)
dist(x, f)

## working with lists
z <- unlist(apply(x, 1, list), recursive = FALSE)
(d <- dist(z))
dist(z, z)

## subsetting
d[[1:2]]
subset(d, c(1,3,4))
d[[c(1,2,2)]]	    # duplicate index gets ignored

## transformations and self-proximities
as.matrix(as.simil(d, function(x) exp(-x)), diag = 1)

## row and column indexes
row.dist(d)
col.dist(d)
}
\keyword{cluster}

