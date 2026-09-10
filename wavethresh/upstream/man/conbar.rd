\name{conbar}
\alias{conbar}
\title{Performs inverse DWT reconstruction step}
\usage{
conbar(c.in, d.in, filter)
}
\arguments{
\item{c.in}{The father wavelet coefficients that you wish to reconstruct in this level's convolution.}
\item{d.in}{The mother wavelet coefficients that you wish to reconstruct in this level's convolution.}
\item{filter}{A given filter that you wish to use in the level reconstruction. This should be the output from the \code{\link{filter.select}} function.}
}
\description{
Wrapper to the C function \code{conbar} which is the main function
in WaveThresh to do filter convolution/reconstruction with data.
Although users use the \code{\link{wr}} function to perform a complete
inverse discrete wavelet transform (DWT) this function repeatedly uses
the \code{conbar} routine, once for each level to reconstruct the next finest
level. The C \code{conbar} routine is possibly the most frequently utilized
by WaveThresh.
}
\details{
The \code{\link{wr}} function performs the inverse wavelet transform on an
\code{\link{wd.object}} class object.

Internally, the \code{\link{wr}} function uses the C \code{conbar} function.
Other functions also make use of \code{conbar} and some R functions also would
benefit from using the fast C code of the \code{conbar} reconstruction hence
this WaveThresh function.
Some of the other functions that use conbar are listed in the SEE ALSO section.
Many other functions call C code that then uses the C version of \code{conbar}.
}
\value{
A vector containing the reconstructed coefficients.
}
\seealso{
\code{\link{av.basis}}
\code{\link{InvBasis.wp}}
\code{\link{wr}}
}
\examples{
#
# Let's generate some test data, just some 32 normal variates.
#
v <- rnorm(32)
#
# Now take the wavelet transform with default filter arguments (which
# are filter.number=10, family="DaubLeAsymm")
#
vwd <- wd(v)
#
# Now, let's take an arbitrary level, say 2, and reconstruct level 3
# scaling function coefficients
#
c.in <- accessC(vwd, lev=2)
d.in <- accessD(vwd, lev=2)
#
conbar(c.in, d.in, filter.select(filter.number=10, family="DaubLeAsymm"))
#[1] -0.50368115  0.04738620 -0.90331807  1.08497622  0.90490528  0.06252717
#[7]  2.55894899 -1.26067508
#
# Ok, this was the pure reconstruction from using only level 2 information.
#
# Let's check this against the "original" level 3 coefficients (which get
# stored on the decomposition step in wd)
#
accessC(vwd, lev=3)
#[1] -0.50368115  0.04738620 -0.90331807  1.08497622  0.90490528  0.06252717
#[7]  2.55894899 -1.26067508
#
# Yep, the same numbers!
#
}
\author{G P Nason}
\keyword{math}
