\name{cl_boot}
\alias{cl_boot}
\title{Bootstrap Resampling of Clustering Algorithms}
\description{
  Generate bootstrap replicates of the results of applying a
  \dQuote{base} clustering algorithm to a given data set.
}
\usage{
cl_boot(x, B, k = NULL,
        algorithm = if (is.null(k)) "hclust" else "kmeans", 
        parameters = list(), resample = FALSE)
}
\arguments{
  \item{x}{the data set of objects to be clustered, as appropriate for
    the base clustering algorithm.}
  \item{B}{an integer giving the number of bootstrap replicates.}
  \item{k}{\code{NULL} (default), or an integer giving the number of
    classes to be used for a partitioning base algorithm.}
  \item{algorithm}{a character string or function specifying the base
    clustering algorithm.}
  \item{parameters}{a named list of additional arguments to be passed to
    the base algorithm.}
  \item{resample}{a logical indicating whether the data should be
    resampled in addition to \dQuote{sampling from the algorithm}.
    If resampling is used, the class memberships of the objects given in
    \code{x} are predicted from the results of running the base
    algorithm on bootstrap samples of \code{x}.}
}
\value{
  A cluster ensemble of length \eqn{B}, with either (if resampling is
  not used, default) the results of running the base algorithm on the
  given data set, or (if resampling is used) the memberships for the
  given data predicted from the results of running the base algorithm on
  bootstrap samples of the data.
}
\details{
  This is a rather simple-minded function with limited applicability,
  and mostly useful for studying the effect of (uncontrolled) random
  initializations of fixed-point partitioning algorithms such as
  \code{\link[stats]{kmeans}} or \code{\link[e1071]{cmeans}}, see the
  examples.  To study the effect of varying control parameters or
  explicitly providing random starting values, the respective cluster
  ensemble has to be generated explicitly (most conveniently by using
  \code{\link{replicate}} to create a list \code{lst} of suitable
  instances of clusterings obtained by the base algorithm, and using
  \code{cl_ensemble(list = lst)} to create the ensemble).
}
\examples{
## Study e.g. the effect of random kmeans() initializations.
data("Cassini")
pens <- cl_boot(Cassini$x, 15, 3)
diss <- cl_dissimilarity(pens)
summary(c(diss))
plot(hclust(diss))
}
\keyword{cluster}
