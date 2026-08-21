\name{rbox}
\alias{rbox}
\title{
Local Regression, Likelihood and Density Estimation.
}
\usage{
rbox(cut=0.8, type="tree", ll, ur)
}
\description{
  \code{rbox()} is used to specify a rectangular box evaluation
  structure for \code{\link{locfit.raw}()}. The structure begins
  by generating a bounding box for the data, then recursively divides
  the box to a desired precision.
}
\arguments{
\item{type}{
If \code{type="tree"}, the cells are recursively divided according to
the bandwidths at each corner of the cell; see Chapter 11 of Loader (1999).
If \code{type="kdtree"}, the K-D tree structure used in Loess
(Cleveland and Grosse, 1991) is used.
}
\item{cut}{
Precision of the tree; a smaller value of \code{cut} results in a
larger tree with more nodes being generated.
}
\item{ll}{
Lower left corner of the initial cell. Length should be the number
of dimensions of the data provided to \code{\link{locfit.raw}()}.
}
\item{ur}{
Upper right corner of the initial cell. By default, \code{ll} and
\code{ur} are generated as the bounding box for the data.
}
}

\examples{
data(ethanol, package="locfit")
plot.eval(locfit(NOx~E+C,data=ethanol,scale=0,ev=rbox(cut=0.8)))
plot.eval(locfit(NOx~E+C,data=ethanol,scale=0,ev=rbox(cut=0.3)))
}

\references{
Loader, C. (1999). Local Regression and Likelihood. Springer, New York.

Cleveland, W. and Grosse, E. (1991). Computational Methods for Local
  Regression. Statistics and Computing 1.
}
\keyword{smooth}
