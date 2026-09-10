\name{image.wd}
\alias{image.wd}
\usage{
\method{image}{wd}(x, strut = 10, type = "D", transform = I, ...)
}
\arguments{
\item{x}{The \code{\link{wd.object}} that you wish to image}
\item{strut}{The width of each coefficient in the image}
\item{type}{Either "C" or "D" depending if you wish to image scaling
	function or wavelet coefficients respectively}
\item{transform}{Apply a numerical transform to the coefficients before
	display}
\item{\dots}{Other arguments}
}
\title{Produce image representation of nondecimated wavelet transform}
\description{
Produces a representation of a nondecimated wavelet transform (time-ordered)
as an image.
}
\details{
Description says all
}
\value{
None
}
\seealso{\code{\link{logabs}}, \code{\link{nlevelsWT}},\code{\link{wd}}}
\examples{
tmp <- wd(rnorm(256), type="station")
\dontrun{image(tmp)}
}
\author{G P Nason}
\keyword{hplot}
