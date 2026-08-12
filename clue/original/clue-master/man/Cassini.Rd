\name{Cassini}
\alias{Cassini}
\title{Cassini Data}
\description{
  A Cassini data set with 1000 points in 2-dimensional space which are
  drawn from the uniform distribution on 3 structures.  The two outer
  structures are banana-shaped; the \dQuote{middle} structure in between
  them is a circle.
}
\usage{data("Cassini")}
\format{
  A classed list with components
  \describe{
    \item{\code{x}}{a matrix with 1000 rows and 2 columns giving the
      coordinates of the points.}
    \item{\code{classes}}{a factor indicating which structure the
      respective points belong to.}
  }
}
\details{
  Instances of Cassini data sets can be created using function
  \code{\link[mlbench]{mlbench.cassini}} in package \pkg{mlbench}.
  The data set at hand was obtained using
  \preformatted{
    library("mlbench")
    set.seed(1234)
    Cassini <- mlbench.cassini(1000)
  }
}
\examples{
data("Cassini")
op <- par(mfcol = c(1, 2))
## Plot the data set:
plot(Cassini$x, col = as.integer(Cassini$classes),
     xlab = "", ylab = "")
## Create a "random" k-means partition of the data:
set.seed(1234)
party <- kmeans(Cassini$x, 3)
## And plot that.
plot(Cassini$x, col = cl_class_ids(party),
     xlab = "", ylab = "")
## (We can see the problem ...)
par(op)
}
\keyword{datasets}
