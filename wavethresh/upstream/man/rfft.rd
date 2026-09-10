\name{rfft}
\alias{rfft}
\title{Real Fast Fourier transform}
\usage{
rfft(x)
}
\arguments{
\item{x}{The vector whose Fourier transform you wish to take}
}
\description{
Compute a real Fast Fourier transform of \code{x}.
}
\details{
Given a vector x this function computes the real continuous Fourier
transform of \code{x}, i.e. it regards \code{x} as points on a periodic
function on [0,1] starting at 0, and finding the coefficients of the functions
1,
\eqn{\sqrt{2}\cos(2\pi t)}{sqrt(2) cos(2 pi t)},
\eqn{\sqrt{2}\sin(2\pi t)}{sqrt(2) sin(2 pi t)}, etc. that gives the expansion
of the interpolant of \code{x}. The number of terms in the expansion 
is the length of \code{x}. If \code{x} is of even length, the last
coefficient will be that of a cosine term with no matching sine.
}
\value{
Returns the Fourier coefficients
}
\seealso{
\code{\link{LocalSpec.wd}},
\code{\link{rfftinv}}}
\examples{
x <- seq(from=0, to=2*pi, length=150)
s1 <- sin(10*x)
s2 <- sin(7*x)
s <- s1 + s2
w <- rfft(s)
\dontrun{ts.plot(w)}
#
# Should see two peaks, corresponding to the two sines at different frequencies
#
}
\author{Bernard Silverman}
\keyword{math}
