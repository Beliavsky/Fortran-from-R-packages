\name{var.add.benchmark}
\alias{var.add.benchmark}
\title{ Expand a Variance Matrix to Include a Benchmark }
\description{
Takes a variance matrix (or 3-dimensional array) and a vector of weights 
for the benchmark.
An object like the input is returned which includes a new asset 
representing the benchmark.
}
\usage{
var.add.benchmark(variance, benchmark.weights, name = "benchmark",
	sum.to.one = TRUE)
}
\arguments{
  \item{variance}{ 
	required.
	A variance matrix or a three-dimensional array 
	where each slice of the third dimension is a variance matrix.
	This must have dimnames which specify the assets, and it must
	include all of the assets named in \code{benchmark.weights}.
	}
  \item{benchmark.weights}{ 
	required.
	A named vector giving the weights of the constituents of the benchmark.
	}
  \item{name}{ 
	character string to be used as the asset name for the benchmark.
	}
  \item{sum.to.one}{ 
	logical vector stating whether to check and enforce that 
	\code{benchmark.weights} sums to one.
	}
}
\value{
a matrix or array similar to the input \code{variance}, but including
an additional asset.
}

\section{Details }{
The (absolute value of) the weights ideally sum to one.
You can give it weights that do not sum to one (perhaps they sum to 100
or the market cap of the benchmark), but you will get a warning that the
weights are being adjusted.

If your "benchmark" is something that doesn't sum to 1 (such as 
portfolio weights minus benchmark weights), then set the \code{sum.to.one}
argument to \code{FALSE}.
}

\section{Revision }{
This help was last revised 2012 January 20.
}

\seealso{ 
\code{\link{var.relative.benchmark}}, \code{\link{var.shrink.eqcor}},
\code{\link{factor.model.stat}}, \code{\link{threeDarr}}.
}

\examples{
\dontrun{
varian.ben <- var.add.benchmark(varian, ftse.constituents, "ftse100")
}

var.orig <- array(c(400, 32, 24, 32, 64, 9.6, 24, 9.6, 144), c(3,3),
	list(c("equities", "bonds", "commodities"), 
	c("equities", "bonds", "commodities")))
var.aa <- var.add.benchmark(var.orig, c(equities=.6, bonds=.4), "e60b40")
}
\keyword{multivariate}
