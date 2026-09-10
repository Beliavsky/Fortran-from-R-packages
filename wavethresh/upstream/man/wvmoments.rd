\name{wvmoments}
\alias{wvmoments}
\title{Compute moments of wavelets or scaling function}
\usage{
wvmoments(filter.number = 10, family = "DaubLeAsymm", moment = 0,
    scaling.function = FALSE)
}
\arguments{
\item{filter.number}{The smoothness of wavelet or scaling function to
	compute moments for, see \code{\link{filter.select}}}
\item{family}{The wavelet family to use, see \code{\link{filter.select}}}
\item{moment}{The moment to compute}
\item{scaling.function}{If \code{FALSE} then a wavelet is used in the
	moment calculation, alternatively if \code{TRUE} the associated
	scaling function is used.}
}

\description{
Numerically compute moments of wavelets or scaling function
}
\details{
Given a wavelet \eqn{\psi(x)}{psi(x)} this function computes the
mth moment \eqn{\int x^m \psi(x) \, dx}{int x^m psi(x) dx}.

Note that for low order moments the integration function often fails for
the usual numerical reasons (this never happened in S!). It might be that
fiddling with the tolerances will improve this situation.
}
\value{
An object of class \code{integrate} containing the integral and other
pieces of interesting information about the moments calculation.
}
\examples{
wvmoments(filter.number=5, family="DaubExPhase", moment=5)
#-1.317600 with absolute error < 7.5e-05
}
\seealso{\code{\link{draw.default}}}
\author{G P Nason}
\keyword{math}
