\name{example.1}
\alias{example.1}
\title{Compute and return piecewise polynomial coordinates.
}
\usage{
example.1()
}
\arguments{
None
}
\description{
This function computes and returns the coordinates of the piecewise polynomial described by Nason and Silverman, 1994. This function is a useful test function for evaluating wavelet shrinkage methodology as it contains smooth parts, a discontinuity and it is periodic.

(Nason, G.P. and Silverman, B.W. (1994) The discrete wavelet transform in
S, \emph{J. Comput. Graph. Statist.}, \bold{3}, 163--191.)
}
\details{
This function computes and returns the x and y coordinates of the piecewise polynomial function described in Nason and Silverman, 1994.
The formula for the piecewise polynomial (which is piecewise cubic) is given in Nason and Silverman, 1994.

The piecewise polynomial returned is a discrete sample on 512 equally spaced points between 0 and 1 (including 0 but excluding 1).

The Donoho and Johnstone test functions can be generated using the
\code{\link{DJ.EX}} function.
}
\value{
A list with two components:
\item{x}{a vector of length 512 containing the ordered x ordinates of the piecewise polynomial.}
\item{y}{a vector of length 512 containing the corresponding y ordinates of the piecewise polynomial.}
}
\seealso{\code{\link{DJ.EX}}}
\examples{
#
# Generate the piecewise polynomial 
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
}
\author{G P Nason}
\keyword{nonparametric}
