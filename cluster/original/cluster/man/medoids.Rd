\name{medoids}
\alias{medoids}
\title{Compute \code{pam}-consistent Medoids from Clustering}
\description{
  Given a data matrix or dissimilarity \code{x} for say \eqn{n}
  observational units and a clustering,
  compute the \code{\link{pam}()}-consistent medoids.
}
\usage{
medoids(x, clustering, diss = inherits(x, "dist"), USE.NAMES = FALSE, ...)
}
\arguments{
  \item{x}{Either a data matrix or data frame, or dissimilarity matrix or
    object, see also \code{\link{pam}}.}
  \item{clustering}{an integer vector of length \eqn{n}, the number of
    observations, giving for each observation the number ('id') of the
    cluster to which it belongs.  In other words, \code{clustering} has
    values from \code{1:k} where \code{k} is the number of clusters, see
    also \code{\link{partition.object}} and \code{\link{cutree}()}, for
    examples where such clustering vectors are computed.}
  \item{diss}{see also \code{\link{pam}}.}
  \item{USE.NAMES}{a logical, typical false, passed to the
    \code{\link{vapply}()} call computing the medoids.}
  \item{\dots}{optional further argument passed to \code{\link{pam}(xj, k=1, \dots)},
    notably \code{metric}, or \code{variant="f_5"} to use a faster algorithm, or
    \code{trace.lev = k}.}
}
%% \details{
%% }
\value{
  a numeric vector of length
}
%% \references{
%% }
\author{Martin Maechler, after being asked how \code{\link{pam}()} could be used
  instead of \code{\link{kmeans}()}, starting from a previous clustering.
}
%% \note{
%% }

\seealso{
  \code{\link{pam}}, \code{\link{kmeans}}.
  Further, \code{\link{cutree}()} and \code{\link{agnes}} (or \code{\link{hclust}}).
}
\examples{
## From example(agnes):
data(votes.repub)
agn1 <- agnes(votes.repub, metric = "manhattan", stand = TRUE)
agn2 <- agnes(daisy(votes.repub), diss = TRUE, method = "complete")
agnS <- agnes(votes.repub, method = "flexible", par.method = 0.625)

for(k in 2:11) {
  print(table(cl.k <- cutree(agnS, k=k)))
  stopifnot(length(cl.k) == nrow(votes.repub), 1 <= cl.k, cl.k <= k, table(cl.k) >= 2)
  m.k <- medoids(votes.repub, cl.k)
  cat("k =", k,"; sort(medoids) = "); dput(sort(m.k), control={})
}

}
\keyword{cluster}
