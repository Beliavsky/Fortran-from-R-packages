\name{lfmarg}
\alias{lfmarg}
\title{
Generate grid margins.
}
\usage{
lfmarg(xlim, m = 40)
}
\arguments{
\item{xlim}{
Vector of limits for the grid. Should be of length 2*d;
the first d components represent the lower left corner,
and the next d components the upper right corner.
Can also be a \code{"locfit"} object.
}
\item{m}{
Number of points for each grid margin. Can be a vector of length d.
}
}
\value{
A list, whose components are the d grid margins.
}
\description{
This function is usually called by \code{\link{plot.locfit}}.
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}
}
\keyword{smooth}
