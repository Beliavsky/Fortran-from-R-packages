\name{cl_tabulate}
\alias{cl_tabulate}
\title{Tabulate Vector Objects}
\description{Tabulate the unique values in vector objects.}
\usage{
cl_tabulate(x)
}
\arguments{
  \item{x}{a vector.}
}
\value{
  A data frame with components:
  \item{values}{the unique values.}
  \item{counts}{an integer vector with the number of times each of the
    unique values occurs in \code{x}.}
}
\examples{
data("Kinship82")
tab <- cl_tabulate(Kinship82)
## The counts:
tab$counts
## The most frequent partition:
tab$values[[which.max(tab$counts)]]
}
\keyword{utilities}
