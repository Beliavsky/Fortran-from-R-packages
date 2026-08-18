\name{loc.min}
\alias{loc.min}
\title{Minimum Location Finder}
\description{This function identifies the location of the minimum value in a square matrix.
}
\usage{
loc.min(my.mat,d)
}
\arguments{
  \item{my.mat}{A square matrix.}
  \item{d}{Dimensions of the matrix.}
}
\value{A vector containing the row and column number of the minimum value.}
\examples{
cmat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
loc.min(my.mat=cmat, d=3)
}
