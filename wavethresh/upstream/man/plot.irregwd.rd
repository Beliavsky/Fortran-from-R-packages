\name{plot.irregwd}
\alias{plot.irregwd}
\title{Plot variance factors of wavelet transform coefficients for irregularly spaced wavelet transform object}
\usage{
\method{plot}{irregwd}(x, xlabels, first.level = 1,
	main = "Wavelet Decomposition Coefficients", scaling = "by.level",
	rhlab = FALSE, sub, ...)
}
\arguments{
\item{x}{The \code{\link{irregwd.objects}}
	object whose coefficients you wish to plot.}
\item{xlabels}{A vector containing the "true" x-axis numbers that went
	with the vector that was transformed to produce the
	irregwd object supplied as the first argument to this function.
	If this argument is missing then the function tries to make up a
	sensible set of x-axis labels.}
\item{first.level}{The first resolution level to begin plotting at.
	This argument can be quite useful when you want to supress some of
	the coarser levels in the diagram.}
\item{main}{The main title of the plot.}
\item{scaling}{How you want the coefficients to be scaled.
	The options are: \code{global} - one scale factor is chosen for
	the whole plot. The scale factor depends on the variance factor to
	be included on the plot that has the largest absolute value.
	The global option is useful when comparing factors that might
	appear anywhere in the plot;
	\code{by.level} - a scale factor is chosen for each resolution level
	in the plot. The scale factor for a level depends on the variance
	factor in that level that has the largest absolute value.
	The \code{by.level} option is useful when you wish to compare
	coefficients within a resolution level.}
\item{rhlab}{If \code{TRUE} then a set of labels is produced on the right
	hand axis. The axis labels in this case refer to the scale factors
	used to scale each level and correspond to value of the largest
	variance factor (in absolute value) in each scale
	(when \code{scaling=="by.level"}) or absolutely (when \code{scaling="global"}). If the \code{rhlab=FALSE} then no right hand axis labels are
	produced.}
\item{sub}{A subtitle for the plot.}
\item{\dots}{Other arguments supplied to the actual plot}
}
\description{
This function plots the variance factors associated with the wavelet coefficients arising from a \code{\link{irregwd.objects}}
irregularly spaced wavelet decomposition object.
}
\details{
Produces a plot similar in style to the ones in Donoho and Johnstone, 1994.
This function is basically the same as
\code{\link{plot.wd}} except that variance factors and not coefficients
are plotted. A variance factor is a number that quantifies the variability of
a coefficient induced by the irregular design that was interpolated to
a regular grid by the \code{\link{makegrid}} function which is used the
by \code{\link{irregwd}} irregular wavelet transform function.

High values of the variance factor correspond to large variance in the
wavelet coefficients but due to the irregular design,
not the original noise structure on the coefficients.
}
\value{
If \code{rhlab==TRUE} then the scaling factors applied to each
scale level are returned. Otherwise \code{NULL} is returned.
}
\examples{
#
# The help for makegrid contains an example
# of using this function.
#
}
\author{Arne Kovac}
\keyword{hplot}
