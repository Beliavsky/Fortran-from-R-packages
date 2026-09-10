\name{plot.wp}
\alias{plot.wp}
\title{Plot wavelet packet transform coefficients}
\usage{
\method{plot}{wp}(x, nvwp = NULL, main = "Wavelet Packet Decomposition", 
    sub, first.level = 5, scaling = "compensated", dotted.turn.on = 5, 
    color.force = FALSE, WaveletColor = 2, NodeVecColor = 3, 
    fast = FALSE, SmoothedLines = TRUE, ...)
}
\arguments{
\item{x}{The wp object whose coefficients you wish to plot.}
\item{nvwp}{An optional associated wavelet packet node vector class object of
	class \code{nvwp}.
	This object is a list of packets in the wavelet packet table.
	If this argument is specified then it is possible to highlight the
	packets in the nvwp objects in a different color using the
	\code{NodeVecColor} argument}
\item{main}{The main title of the plot.}
\item{sub}{A subtitle for the plot.}
\item{first.level}{The first resolution level to begin plotting at.
	This argument can be quite useful when you want to supress some
	of the coarser levels in the diagram.}
\item{scaling}{How you want the coefficients to be scaled.
	The options are: \code{global} - one scale factor is chosen for
	the whole plot. The scale factor depends on the coefficient to be
	included on the plot that has the largest absolute value.
	The global option is useful when comparing coefficients that might
	appear anywhere in the plot; \code{by.level} - a scale factor is
	chosen for each resolution level in the plot.
	The scale factor for a level depends on the coefficient in that level
	that has the largest absolute value. The \code{by.level} option is
	useful when you wish to compare coefficients within a resolution level.
	The other option is \code{compensated} which is the same as global
	except for that finer scales' coefficients are scaled up by a factor
	of SQRT(2)

	I don't know why compensated is the default option? That is probably
	silly!}
\item{dotted.turn.on}{The plot usually includes some dotted vertical bars that
	separate wavelet packets to make it clearer which packets are which.
	This option controls the coarsest resolution level at which dotted
	lines appear. All levels equal to and finer than this level will
	receive the vertical dotted lines.}
\item{color.force}{If FALSE then some "clever" code in CanUseMoreThanOneColor tries to figure out how many colours can be used (THIS HAS NOT BEEN MADE TO WORK
	IN R)  and hence whether colour can be used to pick out wavelet
	packets or elements of a node vector.
	This option was designed to work with S.
	It doesn't work with R and so it is probably best to set
	\code{color.force=T}. In this way no interrogation is done and the
	lines/packets are plotted in the appropriate colours with no
	questions asked.}
\item{WaveletColor}{A colour specification for the colour for wavelet
	coefficients. Wavelet coefficients are a component of wavelet packet
	coefficients and this option allows them to be drawn in a different
	color. In R you can use names like "red", "blue" to select the colors.
	In R you'll also need to set the color.force option to TRUE.}
\item{NodeVecColor}{If a nvwp object is supplied this option can force
	coefficients that are part of that nvwp to be drawn in the specified
	color. See the explanation for the \code{WaveletColor} option above
	about specification in R.}
\item{fast}{This option no longer does anything.}
\item{SmoothedLines}{If TRUE then the scaling function coefficients are
	drawn using lines (and look like mini versions of the original).
	If FALSE then the scaling function coefficients are drawn using
	the \code{segments} function and look like a coarser shadowy
	version of the original.}
\item{\dots}{Other arguments to the plot command}
}
\description{
This function plots wavelet packet transform coefficients arising from a
	\code{\link{wp.object}} object.
}
\details{
A wavelet packet object contains wavelet packet coefficients of a signal
(usually obtained by the \code{\link{wp}} wavelet packet transform function).
Given a wavelet packet object wp it possesses \code{nlevelsWT(wp)}
resolution levels.
In WaveThresh the coarsest level is level 0 and the finest is level
nlevelsWT-1.
For wavelet packets the number of packets at level j is 2^(nlevelsWT-j).

This function plots the wavelet packet coefficients. At the bottom of the plot
the original input function (if present) is plotted.
Then levels above the original plot successively coarser wavelet packet
coefficients. From the Mallat transform point of view smoothing goes up off
the the left of the picture and detail to the right.
The packets are indexed from 0 to the number of packets going from left to
right within each resolution level.

The function has the ability to draw wavelet coefficients in a different color using the \code{WaveletColor} argument.

Optionally, if a node vector wavelet packet object is also supplied, which
contains the specification of a basis selected from the packet table,
then packets in that node vector can be highlighted in a another colour determined by the
\code{NodeVecColor}.

Packets are drawn on the plot and can be separated by vertical dotted lines.
The resolution levels at which this happens can be controlled by the
\code{dotted.turn.on} option.
The coarsest resolution level to be drawn is controlled by the
\code{first.level} option.
}
\value{
Nothing 
}
\seealso{\code{\link{MaNoVe}}, \code{\link{wp}}, \code{\link{wp.object}}}
\examples{
#
# Generate some test data
#
v <- DJ.EX()$blocks
#
# Let's plot these to see what they look like
#
\dontrun{plot(v, type="l")}
#
# Do a wavelet packet transform
#
vwp <- wp(v)
#
# And create a node vector
#
vnv <- MaNoVe(vwp)
#
# Now plot the wavelet packets with the associated node vector
#
\dontrun{plot(vwp, vnv, color.force=T, WaveletColor="red", dotted.turn.on=7)}
#
# The wavelet coefficients are plotted in red. Packets from the node vector
# are depicted in green. The node vector gets plotted after the wavelet
# coefficients so the green packets overlay the red (retry the plot command
# but without the vnv object to see just the
# wavelet coefficients). The vertical dotted lines start at resolution
# level 7.
#
#
}
\author{G P Nason}
\keyword{hplot}
