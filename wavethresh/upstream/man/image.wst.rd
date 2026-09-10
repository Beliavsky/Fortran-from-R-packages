\name{image.wst}
\alias{image.wst}
\usage{
\method{image}{wst}(x, nv, strut = 10, type = "D", transform = I, \dots)
}
\arguments{
\item{x}{The wst object you wish to image}
\item{nv}{An associated node vector, this argument is no longer used
	and should be omitted (in the S version it permitted coloration
	of particular bases)}
\item{strut}{The number of pixels/width that each coefficient should be
	drawn with}
\item{type}{Either "C" or "D" depending on whether you wish to image
	scaling function coefficients or wavelet ones}
\item{transform}{A numerical transform you wish to apply to the coefficients
	before imaging}
\item{\dots}{Other arguments}
}
\title{Produce image representation of a wst class object}
\description{
Produces an image representation of the coefficients contained within
a \code{\link{wst.object}} class object.
}
\details{
Description says all
}
\value{
None
}
\seealso{\code{\link{logabs}},\code{\link{wst}}}
\examples{
tmp <- wst(rnorm(1024))
\dontrun{image(tmp)}
\dontrun{image(tmp, transform=logabs)}
}
\author{G P Nason}
\keyword{hplot}
