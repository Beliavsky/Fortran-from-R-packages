\name{addpkt}
\alias{addpkt}
\title{Add a wavelet packet box to an already set up time-frequency plot}
\usage{
addpkt(level, index, density, col, yvals)
}
\arguments{
\item{level}{The level at which the box or yvals are plotted}
\item{index}{The packet index at which the box of yvals are plotted}
\item{density}{The density of the shading of the box}
\item{col}{The color of the box}
\item{yvals}{If this argument is missing then a shaded coloured box is
	drawn, otherwise a time series of \code{yvals} is plotted where
	the box would have been.}
}
\description{
This function assumes that a high-level plot has already been set
up using \code{\link{plotpkt}}. Given that this function plots
a wavelet packet box at a given level, packet index and with
particular shading and color and optionally plotting a sequence of
coefficients at that location rather than a shaded box.
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{basisplot}},\code{\link{basisplot.BP}},
	\code{\link{basisplot.wp}}, \code{\link{plotpkt}}, \code{\link{plot.nvwp}}}
\author{G P Nason}
\keyword{dplot}
