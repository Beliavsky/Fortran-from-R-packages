\name{gsiplotmargin}
\alias{gsi.plotmargin}
\title{Internal functions of the compositions package}
\description{
  Internal function to compute 2D marginal compositions for plots  
}
\usage{
gsi.plotmargin(X,d,margin,what="data")
}
\arguments{
  \item{X}{a multivariate compositional dataset}
  \item{d}{a numeric or character vector of two elements specifying the
    parts to be kept in the marginalization process}
  \item{margin}{a character specifying the type of margin to be
    chosen. Possible values are "acomp", "rcomp" or a column name from
    the dataset.}
  \item{what}{simply here to avoid an error if the argument is passed.}
}
\value{
  a composition of three elements.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{gsi}}
}
\keyword{internal}
