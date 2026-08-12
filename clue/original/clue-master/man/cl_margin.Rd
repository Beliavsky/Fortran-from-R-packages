\name{cl_margin}
\alias{cl_margin}
\title{Membership Margins}
\description{
  Compute the \emph{margin} of the memberships of a partition, i.e., the
  difference between the largest and second largest membership values of
  the respective objects.
}
\usage{
cl_margin(x)
}
\arguments{
  \item{x}{an \R object representing a partition of objects.}
}
\details{
  For hard partitions, the margins are always 1.

  For soft partitions, the margins may be taken as an indication of the
  \dQuote{sureness} of classifying an object to the class with maximum
  membership value.
}
\examples{
data("GVME")
## Look at the classes obtained for 1980:
split(cl_object_names(GVME[["1980"]]), cl_class_ids(GVME[["1980"]]))
## Margins:
x <- cl_margin(GVME[["1980"]])
## Add names, and sort:
names(x) <- cl_object_names(GVME[["1980"]])
sort(x)
## Note the "uncertainty" of assigning Egypt to the "intermediate" class
## of nations.
}
\keyword{cluster}
