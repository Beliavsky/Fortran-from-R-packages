\name{threeDarr}
\alias{threeDarr}
\title{ 
	Combine matrices into 3D array
}
\description{
	Returns a three-dimensional array given one or more matrices
	and instructions on how to combine them.
}
\usage{
threeDarr(..., rep = 1, union = TRUE, slicenames = NULL)
}
\arguments{
  \item{\dots}{ 
	one or more matrices.
}
  \item{rep}{ 
	an integer saying how many times to replicate the slices.
}
  \item{union}{ 
	logical value: if \code{TRUE}, then the union of row and column
	names is created in the output.
	If \code{FALSE}, then the intersection is done.
}
  \item{slicenames}{ 
	a character vector with length equal to the number of slices
	in the resulting array.
}
}
\value{
	an array with the first two dimensions being the union or
	intersection of the first two dimensions of the input matrices,
	and third dimension equal to the number of input matrices times
	the number of replications.
}

\section{Details }{
The full name of arguments \code{rep}, \code{union} and \code{slicenames}
must be given (no abbreviations)
because they come after the three-dots construct.
}

\section{Revision }{
This help was last revised 2012 January 22.
}

\seealso{ 
	More general functionality of this sort can be found in the
	\code{abind} package.
}
\examples{
\dontrun{
multiple.var <- threeDarr(var1, var2, slicenames=c("standardVar", "crashVar"))

reparr <- threeDarr(matrix1, rep=3)
}
}
\keyword{array }
