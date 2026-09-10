\name{sure}
\alias{sure}
\usage{
sure(x)
}
\arguments{
\item{x}{Vector of (normalized) wavelet coefficients. Coefficients
	should be supplied divided by their standard deviation, or
	some robust measure of scale}
}
\title{Computes the minimum of the SURE thresholding function}
\description{
Computes the minimum of the SURE thresholding function for
wavelet shrinkage as described in Donoho, D.L. and Johnstone, I.M.
(1995) Adapting to unknown smoothness via wavelet shrinkage.
\emph{J. Am. Statist. Ass.}, \bold{90}, 1200-1224.
}
\details{
SURE is a method for unbiasedly estimating the risk of an estimator. Stein (1981) showed that for a nearly arbitrary, nonlinear biased estimator, one can estimate its loss unbiasedly. See the Donoho and Johnstone, 1995 for further references and explanation. This function minimizes formula (11) from that paper.
}
\seealso{\code{\link{threshold}}}
\value{
The absolute value of the wavelet coefficient that minimizes the
SURE criteria 
}
\examples{
#
# Let's create "pretend" vector of wavelet coefficients contaminated with
# "noise".
#
v <- c(0.1, -0.2, 0.3, -0.4, 0.5, 99, 12, 6)
#
# Now, what's sure of this?
#
sure(v)
#
# [1] 0.5
#
#
# I.e. the large significant coefficients are 99, 12, 6 and the noise is
# anything less than this in abs value. So sure(v) is a good point to threshold
# at.
}
\author{G P Nason}
\keyword{math}
