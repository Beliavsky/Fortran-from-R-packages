\name{cl_medoid}
\alias{cl_medoid}
\title{Medoid Partitions and Hierarchies}
\description{
  Compute the medoid of an ensemble of partitions or hierarchies, i.e.,
  the element of the ensemble minimizing the sum of dissimilarities to
  all other elements.
}
\usage{
cl_medoid(x, method = "euclidean")
}
\arguments{
  \item{x}{an ensemble of partitions or hierarchies, or something
    coercible to that (see \code{\link{cl_ensemble}}).}
  \item{method}{a character string or a function, as for argument
    \code{method} of function \code{\link{cl_dissimilarity}}.}
}
\value{
  The medoid partition or hierarchy.
}
\details{
  Medoid clusterings are special cases of \dQuote{consensus} clusterings
  characterized as the solutions of an optimization problem.
  See \bibcitet{R:Gordon:1999} for more information.
  
  The dissimilarities \code{d} for determining the medoid are obtained
  by calling \code{cl_dissimilarity} with arguments \code{x} and
  \code{method}.  The medoid can then be found as the (first) row index
  for which the row sum of \code{as.matrix(d)} is minimal.  Modulo
  possible differences in the case of ties, this gives the same results
  as (the medoid obtained by) \code{\link[cluster]{pam}} in package
  \pkg{cluster}.
}
\references{
  \bibshow{R:Gordon:1999}
}  
\seealso{
  \code{\link{cl_consensus}}
}
\examples{
## An ensemble of partitions.
data("CKME")
pens <- CKME[1 : 20]
m1 <- cl_medoid(pens)
diss <- cl_dissimilarity(pens)
require("cluster")
m2 <- pens[[pam(diss, 1)$medoids]]
## Agreement of medoid consensus partitions.
cl_agreement(m1, m2)
## Or, more straightforwardly:
table(cl_class_ids(m1), cl_class_ids(m2))
}
\keyword{cluster}
