\name{plotdenwd}
\alias{plotdenwd}
\title{Plot the wavelet coefficients of a p.d.f.}
\usage{
plotdenwd(wd, xlabvals, xlabchars, ylabchars, first.level=0,
	top.level=nlevelsWT(wd)-1,
	main="Wavelet Decomposition Coefficients", scaling="global",
	rhlab=FALSE, sub, NotPlotVal=0.005, xlab="Translate",
	ylab="Resolution Level", aspect="Identity", \dots)
}
\arguments{
\item{wd}{Wavelet decomposition object, usually output from \code{\link{denwd}},
	possibly thresholded.}
\item{xlabvals}{X-axis values at which the \code{xlabchars} will be printed}
\item{xlabchars}{The x-label characters to be plotted at \code{xlabvals}}
\item{ylabchars}{The y-label characters}
\item{first.level}{This specifies how many of the coarse levels of
	coefficients are omitted from the plot.
	The default value of 0 means that all levels are plotted.}
\item{top.level}{This tells the plotting rountine the true resolution level of
	the finest level of coefficients. The default results in the
	coarsest level being labelled 0. The "correct" value can be
	determined from the empirical scaling function coefficient object
	(output from denproj) as in the example below.}
\item{main}{The title of the plot.}
\item{scaling}{The type of scaling applied to levels within the plot.
	This can be "compensated", "by.level" or "global".
	See \code{\link{plot.wd}} for further details.}
\item{rhlab}{Determines whether the scale factors applied to each level
	before plotting are printed as the right hand axis.}
\item{sub}{The plot subtitle}
\item{NotPlotVal}{If the maximum coefficient in a particular level is smaller than \code{NotPlotVal}, then the level is not plotted.}
\item{xlab}{The x-axis label}
\item{ylab}{The y-axis label}
\item{aspect}{Function to apply to coefficients before plotting}
\item{\dots}{Other arguments to the main plot routine}
}

\description{
Plots the wavelet coefficients of a density function.
}
\details{
Basically the same as
\code{\link{plot.wd}} except that it copes with the zero boundary conditions
used in density estimation. Note that for large filter number wavelets the
high level coefficients will appear very squashed compared with the low
level coefficients. This is a consequence of the zero boundary conditions
and the use of the convention that each coefficient is plotted midway between
two coefficients at the next highest level, as in \code{\link{plot.wd}}.
}
\value{
Axis labels to the right of the picture (scale factors).
These are returned as they are sometimes hard to read on the plot.
}
\examples{
# Simulate data from the claw density, find the empirical
# scaling function coefficients, decompose them and plot
# the resulting wavelet coefficients

data <- rclaw(100)
datahr <- denproj(data, J=8, filter.number=2, family="DaubExPhase")
data.wd <- denwd(datahr)
\dontrun{plotdenwd(data.wd, top.level=(datahr$res$J-1))}
#
# Now use a smoother wavelet
#
datahr <- denproj(data, J=8, filter.number=10, family="DaubLeAsymm")
data.wd <- denwd(datahr)
\dontrun{plotdenwd(data.wd, top.level=(datahr$res$J-1))}
}
\author{David Herrick}
\keyword{hplot}
