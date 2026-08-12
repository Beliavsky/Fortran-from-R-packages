\name{cl_membership}
\alias{cl_membership}
\alias{as.cl_membership}
\title{Memberships of Partitions}
\description{
  Compute the memberships values for objects representing partitions.
}
\usage{
cl_membership(x, k = n_of_classes(x))
as.cl_membership(x)
}
\arguments{
  \item{x}{an R object representing a partition of objects (for
    \code{cl_membership}) or raw memberships or class ids (for
    \code{as.cl_membership}).}
  \item{k}{an integer giving the number of columns (corresponding to
    class ids) to be used in the membership matrix.  Must not be less,
    and default to, the number of classes in the partition.}
}
\value{
  An object of class \code{"cl_membership"} with the matrix of
  membership values.
}
\details{
  \code{cl_membership} is a generic function.

  The methods provided in package \pkg{clue} handle the partitions
  obtained from clustering functions in the base R distribution, as well
  as packages \pkg{RWeka}, \pkg{cba}, \pkg{cclust}, \pkg{cluster},
  \pkg{e1071}, \pkg{flexclust}, \pkg{flexmix}, \pkg{kernlab},
  \pkg{mclust}, \pkg{movMF} and \pkg{skmeans} (and of course, \pkg{clue}
  itself).

  \code{as.cl_membership} can be used for coercing \dQuote{raw} class
  ids (given as atomic vectors) or membership values (given as numeric
  matrices) to membership objects.
}
\seealso{
  \code{\link{is.cl_partition}}
}
\examples{
## Getting the memberships of a single soft partition.
d <- dist(USArrests)
hclust_methods <-
    c("ward", "single", "complete", "average", "mcquitty")
hclust_results <- lapply(hclust_methods, function(m) hclust(d, m))
names(hclust_results) <- hclust_methods 
## Now create an ensemble from the results.
hens <- cl_ensemble(list = hclust_results)
## And add the results of agnes and diana.
require("cluster")
hens <- c(hens, list(agnes = agnes(d), diana = diana(d)))
## Create a dissimilarity object from this.
d1 <- cl_dissimilarity(hens)
## And compute a soft partition.
party <- fanny(d1, 2)
round(cl_membership(party), 5)
## The "nearest" hard partition to this:
as.cl_hard_partition(party)
## (which has the same class ids as cl_class_ids(party)).

## Extracting the memberships from the elements of an ensemble of
## partitions.
pens <- cl_boot(USArrests, 30, 3)
pens
mems <- lapply(pens, cl_membership)
## And turning these raw memberships into an ensemble of partitions.
pens <- cl_ensemble(list = lapply(mems, as.cl_partition))
pens
pens[[length(pens)]]
}
\keyword{cluster}
