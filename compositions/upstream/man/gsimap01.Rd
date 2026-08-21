\name{gsi.mapin01}
\alias{gsi.mapin01}
\alias{gsi.mapfrom01}
\alias{gsi.mapmin}
\alias{gsi.mapmax}
\title{Internal functions: Storing integers as reals}
\description{
  An integer number is stored in a dataset with a given range.
}
\usage{
gsi.mapin01(i,min=0,max=1)
gsi.mapfrom01(x)
gsi.mapmin(x)
gsi.mapmax(x)
}
\arguments{
  \item{i}{an integer number to be masked}
  \item{x}{a numeric vector created by \code{gsi.mapin01}.}
  \item{max}{the maximum of the created dataset}
  \item{min}{the minimum of the created dataset}
}
\value{
  \code{gsi.mapin01} gives a vector \var{x} with
  \code{range(\var{x})==c(min,max)} and \code{gsi.mapfrom01(\var{x})},
  \code{gsi.mapmin(\var{x})}, \code{gsi.mapmax(\var{x})} result in \code{i},
  \code{max} and \code{min}. 
}
\details{
  The function is used to get full control over the graphic ranges in
  pair plots and to pass the used column to panel functions.
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{gsi}}, \code{\link{plot.acomp}}
}
\examples{
#gsi.mapin01(5)
}
\keyword{internal}
