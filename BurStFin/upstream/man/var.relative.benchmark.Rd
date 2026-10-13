\name{var.relative.benchmark}
\alias{var.relative.benchmark}
\title{ Transform a Variance Matrix to be Relative to a Benchmark }
\description{
Returns a matrix (or 3-dimensional array) with dimensions 
one smaller than the input.
The returned variance is relative to the benchmark.
}
\usage{
var.relative.benchmark(variance, benchmark)
}
\arguments{
  \item{variance}{ 
	required.
	A variance matrix or a three-dimensional array 
	where each slice of the third dimension is a variance matrix.
	}
  \item{benchmark}{ 
	required.
	A character string naming which asset in \code{variance} is to be
	used as the benchmark.
	}
}
\value{
a matrix or array similar to the input \code{variance}, but with one
less asset and containing variances that are relative to that asset.

There is a \code{call} attribute which gives the command that created
the object (and hence the benchmark to which the object is relative).
}

\section{Revision }{
This help was last revised 2012 January 20.
}

\seealso{ 
\code{\link{var.add.benchmark}}, \code{\link{var.shrink.eqcor}},
\code{\link{factor.model.stat}}, \code{\link{threeDarr}}.
}

\examples{
var.orig <- array(c(400, 32, 24, 32, 64, 9.6, 24, 9.6, 144), c(3,3),
	list(c("equities", "bonds", "commodities"), 
	c("equities", "bonds", "commodities")))
var.aa <- var.add.benchmark(var.orig, c(equities=.6, bonds=.4), "e60b40")

var.rel <- var.relative.benchmark(var.aa, "e60b40")
}
\keyword{multivariate}
