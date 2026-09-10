\name{plot.wst}
\alias{plot.wst}
\title{Plot packet-ordered non-decimated wavelet transform coefficients.
}
\usage{
\method{plot}{wst}(x, main = "Nondecimated Wavelet (Packet) Decomposition", 
    sub, first.level = 5, scaling = "compensated", dotted.turn.on = 5, 
    aspect = "Identity", ...) 
}
\arguments{
\item{x}{The wst object whose coefficients you wish to plot.}
\item{main}{The main title of the plot.}
\item{sub}{A subtitle for the plot.}
\item{first.level}{The first resolution level to begin plotting at.
	This argument can be quite useful when you want to suppress some of
	the coarser levels in the diagram.}
\item{scaling}{How you want the coefficients to be scaled.
	The options are: \code{global} - one scale factor is chosen for the
	whole plot. The scale factor depends on the coefficient to be
	included on the plot that has the largest absolute value.
	The global option is useful when comparing coefficients that might
	appear anywhere in the plot;
	\code{by.level} - a scale factor is chosen for each resolution level
	in the plot. The scale factor for a level depends on the coefficient
	in that level that has the largest absolute value. The
	\code{by.level} option is useful when you wish to compare coefficients
	within a resolution level.
	The other option is \code{compensated} which is the same as global
	except for that finer scales' coefficients are scaled up by a factor
	of SQRT(2)

	I don't know why compensated is the default option? It is a bit silly.}
\item{dotted.turn.on}{The plot usually includes some dotted vertical bars that
	separate wavelet packets to make it clearer which packets are which.
	This option controls the coarsest resolution level at which dotted
	lines appear. All levels equal to and finer than this level will
	receive the vertical dotted lines.}
\item{aspect}{A transform to apply to the coefficients before plotting.
	If the coefficients are complex-valued and aspect="Identity" then
	the modulus of the coefficients are plotted.}
\item{\dots}{Other arguments to plot}
}
\description{
This function plots packet-ordered non-decimated wavelet transform
coefficients arising from a \code{\link{wst.object}} object.
}
\details{
A packet-ordered non-decimated wavelet object contains coefficients
of a signal (usually obtained by the \code{\link{wst}}
packet-ordered non-decimated wavelet transform, but also
functions that derive such objects, such as \code{\link{threshold.wst}}).

A packet-ordered nondecimated wavelet object, x,
possesses \code{nlevelsWT(x)} resolution levels.
In WaveThresh the coarsest level is level 0 and the finest is level
\code{nlevelsWT-1}. For packet-ordered nondecimated wavelet
the number of blocks (packets) at
level \code{j} is \code{2^(nlevelsWT-j)}.

This function plots the coefficients.
At the bottom of the plot the original input function (if present) is plotted.
Then levels above the original plot successively coarser wavelet 
coefficients.
Each packet of coefficients is plotted within dotted vertical lines.
At the finest level there are two packets: one (the left one) correspond to
the wavelet coefficients that would be obtained using the (standard)
decimated wavelet transform function, \code{\link{wd}},
and the other packet are those coefficients that would have been obtained
using the standard decimated wavelet transform after a unit cyclic shift.

For coarser levels there are more packets corresponding to different cyclic
shifts (although the computation is not performed using shifting operations
the effect is the same). For full details see Nason and Silverman, 1995.

Packets are drawn on the plot and can be separated by vertical dotted lines.
The resolution levels at which this happens can be controlled by the
\code{dotted.turn.on} option.
The coarsest resolution level to be drawn is controlled by
the \code{first.level option}.

\emph{It should be noted that the packets referred to here are just the
	blocks of nondecimated wavelet coefficients in a packet-ordering.
	These are different to wavelet packets (produced by \code{\link{wp}})
	and nondecimated wavelet packets (produced by \code{\link{wpst}})}
}
\value{
Nothing
}
\seealso{
\code{\link{MaNoVe}},\code{\link{threshold.wst}}, \code{\link{wst}}, \code{\link{wst.object}}}
\examples{
#
# Generate some test data
#
v <- DJ.EX()$heavi
#
# Let's plot these to see what they look like
#
\dontrun{plot(v, type="l")}
#
# Do a packet-ordered non-decimated wavelet packet transform
#
vwst <- wst(v)
#
# Now plot the coefficients
#
\dontrun{plot(vwst)}
#
# Note that the "original" function is at the bottom of the plot.
# The finest scale coefficients (two packets) are immediately above.
# Increasingly coarser scale coefficients are above that!
#
}
\author{G P Nason}
\keyword{hplot}
