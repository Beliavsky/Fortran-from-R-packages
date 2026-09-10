\name{plot.imwd}
\alias{plot.imwd}
\alias{plot.imwdc}
\title{Draw a picture of the 2D wavelet coefficients using image}
\usage{
\method{plot}{imwd}(x, scaling = "by.level", co.type = "abs",
	package = "R", plot.type = "mallat", arrangement = c(3, 3),
	transform = FALSE, tfunction = sqrt, ...)
\method{plot}{imwdc}(x, verbose=FALSE, ...)
}
\arguments{
\item{x}{The 2D imwd object you wish to depict}
\item{scaling}{How coefficient scaling is performed. The options
	are \code{by.level} to scale the coefficients independently
	by level, anything else causes coefficients to be scaled
	globally}
\item{co.type}{Can be \code{"abs"} for the absolute values of the
	coefficients to be plotted, can be \code{"mabs"} for the
	negative absolute values or \code{"none"} for none of this.}
\item{package}{Can be \code{"R"} for the R package, or \code{"S"}.
	The latter does less interesting things and results in a simpler
	plot}
\item{plot.type}{If this argument is \code{"mallat"} the coefficients
	at different scales and orientations are packed into one image
	and plotted, a format originating from Mallat's early papers on this.
	The other possibility is \code{"cols"} which plots each combination
	of scale and direction on a separate plot. This latter format is
	useful for examining coefficients, especially at the coarser scales.}
\item{arrangement}{If \code{plot.type="cols"} then this argument
	specifies how many rows and columns there are in the plot
	array.}
\item{transform}{If FALSE then the coefficients are plotted as they
	are (subject to the \code{co.type} argument above), if TRUE
	then the transform function supplied by \code{tfunction} is
	applied to the coefficients.}
\item{tfunction}{If \code{transform=TRUE} then this function gets applied
	to transform the coefficients before plotting}
\item{verbose}{Print out informative messages}
\item{...}{Supply other arguments to the call to the \code{image} function.
	This is very useful to, e.g., can the colours, or other aspects of
	the image}
}

\description{
This function images 2D the absolute values discrete wavelet transform coefficients arising from a \code{\link{imwd.object}} object.
}
\details{
Description says all
}
\value{
If the \code{package="S"} argument is set then a matrix is returned
containing the image that would have been plotted (and this only works
if the \code{plot.type="mallat"} argument is set also.
}
\seealso{\code{\link{imwd}}, \code{\link{imwd.object}}, \code{\link{threshold.imwd}}}
\examples{
data(lennon)
lwd <- imwd(lennon)
\dontrun{plot(lwd)}
\dontrun{plot(lwd, col=grey(seq(from=0, to=1, length=100)), transform=TRUE)}
}
\author{G P Nason}
\keyword{hplot}
