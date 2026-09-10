\name{MaNoVe.wp}
\alias{MaNoVe.wp}
\title{Make Node Vector (using Coifman-Wickerhauser best-basis type algorithm) on wavelet packet object}
\usage{
\method{MaNoVe}{wp}(wp, verbose=FALSE, \dots)
}
\arguments{
\item{wp}{The wp object for which you wish to find the best basis for.}
\item{verbose}{Whether or not to print out informative messages}
\item{\dots}{Other arguments}
}
\description{
This method chooses a "best-basis" using the Coifman-Wickerhauser (1992)
algorithm applied to wavelet packet, \code{\link{wp.object}}, objects.
}
\details{
Description says all
}
\value{
A wavelet packet node vector object of class \code{nvwp},
a basis description. This can
be fed into a basis inversion using, say, the function
\code{\link{InvBasis}}.
}
\seealso{
\code{\link{InvBasis}},
\code{\link{MaNoVe}},
\code{\link{MaNoVe.wst}},
\code{\link{wp.object}},
\code{\link{wp}}
}
\examples{
#
# See example of use of this function in the examples section
# of the help of plot.wp
#
# A node vector vnv is created there that gets plotted.
#
}
\author{G P Nason}
\keyword{smooth}
