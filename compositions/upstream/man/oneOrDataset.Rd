\name{oneOrDataset}
\alias{oneOrDataset}
\title{Treating single compositions as one-row datasets}
\description{
  A dataset is converted to a data matrix. A single data item (i.e. a
  simple vector) is converted to a one-row data matrix. 
}
\usage{
          oneOrDataset(W,B=NULL)
          }
\arguments{
  \item{W}{a vector, matrix or dataframe}
  \item{B}{an optional second vector, matrix or data frame having the
    intended number of rows.}
}
\value{
   A data matrix containing the same data as W. If W is a vector it is
   interpreded as a single row. If \code{B} is given and
   \code{length(dim(B))!= 2} and \code{W} is a vector,
   then \code{W} is repeated \code{nrow(B)} times.
 }
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
oneOrDataset(c(1,2,3))
oneOrDataset(c(1,2,3),matrix(1:12,nrow=4))
oneOrDataset(data.frame(matrix(1:12,nrow=4)))
}
\keyword{multivariate}
