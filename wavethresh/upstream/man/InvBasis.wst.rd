\name{InvBasis.wst}
\alias{InvBasis.wst}
\title{Invert a wst library representation with a basis specification}
\usage{
\method{InvBasis}{wst}(wst, nv, \dots)
}
\arguments{
\item{wst}{The wst object that you wish to invert}
\item{nv}{The node vector, basis spec, that you want to pick out}
\item{...}{Other arguments, that don't do anything here}
}
\description{
Inverts a wst basis representation with a given basis specification,
for example an output from the \code{\link{MaNoVe}} function.
}
\details{
Objects arising from a \code{\link{wst.object}} specification
are a representation of a signal with respect to a library
of basis functions. A particular basis specification can be obtained
using the \code{\link{numtonv}} function which can pick an indexed
basis function, or \code{\link{MaNoVe.wst}} which uses the
Coifman-Wickerhauser minimum entropy method to select a basis.
This function takes a \code{\link{wst.object}} and
a particular basis description (in a \code{\link{nv.object}} node vector
object) and inverts the representation with respect to that selected basis.
}
\value{
The inverted reconstruction
}
\seealso{\code{\link{numtonv}},\code{\link{nv.object}},\code{\link{MaNoVe.wst}},\code{\link{threshold.wst}},\code{\link{wst}}}
\examples{
#
# Let's generate a noisy signal
#
x <- example.1()$y + rnorm(512, sd=0.2)
#
# You can plot this if you like
#
\dontrun{ts.plot(x)}
#
# Now take the nondecimated wavelet transform
#
xwst <- wst(x)
#
# Threshold it
#
xwstT <- threshold(xwst)
#
# You can plot this too if you like
#
\dontrun{plot(xwstT)}
#
# Now use Coifman-Wickerhauser to get a "good" basis
#
xwstTNV <- MaNoVe(xwstT)
#
# Now invert the thresholded wst using this basis specification
#
xTwr <- InvBasis(xwstT, xwstTNV)
#
# And plot the result, and superimpose the truth in dotted
#
\dontrun{ts.plot(xTwr)}
\dontrun{lines(example.1()$y, lty=2)}
}
\author{G P Nason}
\keyword{smooth}
