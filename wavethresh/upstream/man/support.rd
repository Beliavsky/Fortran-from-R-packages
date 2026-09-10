\name{support}
\alias{support}
\title{Returns support of compactly supported wavelets.}
\usage{
support(filter.number=10, family="DaubLeAsymm", m=0, n=0)
}
\arguments{
\item{filter.number}{The member index of a wavelet within the family.
	For Daubechies' compactly supported wavelet this is the number of
	vanishing moments which is related to the smoothness.
	See \code{\link{filter.select}} for more information on the
	wavelets.}
\item{family}{The family of wavelets.
	See \code{\link{filter.select}} for more information on the wavelets.}
\item{m}{Optional scale value (in usual wavelet terminology this is j)}
\item{n}{Optional translation value (in usual wavelet terminology, this is k)}
}

\description{
Returns the support for compactly supported wavelets. This information is useful for drawing wavelets for annotating axes.
}
\details{
It is useful to know the support of a wavelet when drawing it to annotate
labels. Other functions, such as wavelet density estimation
(\code{\link{CWavDE}}), also use this information.
}
\value{
A list with the following components (each one is a single numeric value)
\item{lh}{Left hand support of the wavelet with scale m and translation n.
	These values change as m and n (although when m=0 the function
	confusingly returns the next coarser wavelet where you might
	expect it to return the mother. The mother is indexed by m=-1)}
\item{rh}{As lh but returns the rh end.}
\item{psi.lh}{left hand end of the support interval for the mother wavelet
	(remains unchanged no matter what m or n are)}
\item{psi.rh}{right hand end of the support interval for the mother wavelet
	(remains unchanged no matter what m or n are)}
\item{phi.lh}{left hand end of the support interval for the father wavelet
	(remains unchanged no matter what m or n are)}
\item{phi.rh}{right hand end of the support interval for the father wavelet (remains unchanged no matter what m or n are)}
}
\seealso{
\code{\link{CWavDE}},
\code{\link{draw.default}},
\code{\link{filter.select}}}
\examples{
#
# What is the support of a Haar wavelet?
#
support(filter.number=1, family="DaubExPhase", m=0, n=0)
#$lh
#[1] 0
#
#$rh
#[1] 2
#
#$psi.lh
#[1] 0
#
#$psi.rh
#[1] 1
#
#$phi.lh
#[1] 0
#
#$phi.rh
#[1] 1
#
# So the mother and father wavelet have support [0,1]
#
}
\author{G P Nason}
\keyword{math}
