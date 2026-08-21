\name{plotbyfactor}
\alias{plotbyfactor}
\title{
x-y scatterplot, colored by levels of a factor.
}
\usage{
plotbyfactor(x, y, f, data, col = 1:10, pch = "O", add = FALSE, lg,
    xlab = deparse(substitute(x)), ylab = deparse(substitute(y)),
    log = "", ...)
}
\description{
  Produces a scatter plot of x-y data, with different classes given
  by a factor f. The different classes are identified by different
  colours and/or symbols.
}
\arguments{
\item{x}{
Variable for x axis.
}
\item{y}{
Variable for y axis.
}
\item{f}{
Factor (or variable for which as.factor() works).
}
\item{data}{
data frame for variables x, y, f. Default: sys.parent().
}
\item{col}{
Color numbers to use in plot. Will be replicated if shorter than the
number of levels of the factor f. Default: 1:10.
}
\item{pch}{
Vector of plot characters. Replicated if necessary.
Default: "O".
}
\item{add}{
If \code{TRUE}, add to existing plot. Otherwise, create new plot.
}
\item{lg}{
Coordinates to place a legend. Default: Missing (no legend).
}
\item{xlab, ylab}{Axes labels.}
\item{log}{Should the axes be in log scale? Use \code{"x"}, \code{"y"},
  or \code{"xy"} to specify which axis to be in log scale.}
\item{...}{
Other graphical parameters, labels, titles e.t.c.
}}
\examples{
data(iris)
plotbyfactor(petal.wid, petal.len, species, data=iris)
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a3.
