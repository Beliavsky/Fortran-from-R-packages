\name{draw.default}
\alias{draw.default}
\title{Draw picture of a wavelet or scaling function. }
\description{
This function can produce pictures of one- or two-dimensional wavelets or scaling functions at various levels of resolution. 
}
\usage{
\method{draw}{default}(filter.number = 10, family = "DaubLeAsymm", resolution = 8192,
    verbose = FALSE, plot.it = TRUE, main = "Wavelet Picture", sub = zwd$
    filter$name, xlab = "x", ylab = "psi", dimension = 1, twodplot
     = persp, enhance = TRUE, efactor = 0.05, scaling.function = FALSE, 
	type="l",
    \dots)
}
\arguments{
\item{filter.number }{This selects the index number of the wavelet or scaling function you want to draw (from within the wavelet family).}
\item{family}{specifies the family of wavelets that you want to draw. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{resolution }{specifies the resolution that the wavelet or scaling function is computed to. It does not necessarily mean that you see all of these points as if the enhance option is TRUE then some function points are omitted.}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
\item{plot.it }{If TRUE then this function attempts to plot the function (i.e. draw it on a graphics device which should be active). If FALSE then this function returns the coordinates of the object that would have been plotted.}
\item{main}{a main title for the plot}
\item{sub}{a subtitle for the plot.}
\item{xlab}{a label string for the x-axis}
\item{ylab}{a label string for the y-axis}
\item{dimension }{whether to make a picture of the one-dimensional wavelet or the two-dimensional wavelet.}
\item{twodplot }{which function to use to produce the two-dimensional plot if dimension=2. The function you supply should accept data just like the contour or persp functions supplied with S-Plus.}
\item{enhance}{If this argument is TRUE then the plot is enhanced in the following way. Many of Daubechies' compactly supported wavelets are near to zero on a reasonable proportion of their support. So if such a wavelet is plotted quite a lot of it looks to be near zero and the interesting detail seems quite small. This function chooses a nice range on which to plot the central parts of the function and the function is plotted on this range.}
\item{efactor}{Variable which controls the range of plotting associated with the enhance option above. Any observations smaller than efactor times the range of the absolute function values are deemed to be too small. Then the largest range of ``non-small'' values is selected to be plotted.}
\item{scaling.function }{If this argument is TRUE the scaling function is plotted otherwise the mother wavelet is plotted.}
\item{type}{The \code{type} argument supplied to the plot command}
\item{\dots}{other arguments you can supply to the plot routine embedded within this function.}

}
\details{
The algorithm underlying this function produces a picture of the wavelet or scaling function by first building a \code{wavelet decomposition} object of the correct size (i.e. \code{correct resolution}) and setting all entries equal to zero. Then one coefficient at a carefully selected resolution level is set to one and the decomposition is inverted to obtain a picture of the wavelet. 
}
\value{
If \code{plot.it=FALSE} then usually a list containing coordinates of the object that \emph{would} have been plotted is returned. This can be useful if you don't want S-Plus to do the plotting or you wish to use the coordinates of the wavelets for other purposes.} 

\note{A plot is produced of the wavelet or scaling function if \code{plot.it=TRUE}.} 
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{filter.select}}, \code{\link{ScalingFunction}},\code{\link{wd}}, \code{\link{wd.object}}, \code{\link{wr}}, \code{\link{wr.wd}}. 
}
\examples{
#
# First make sure that your favourite graphics device is open
# otherwise S-Plus will complain.
#
# Let's draw a one-dimensional Daubechies least-asymmetric wavelet
# N=10
#
\dontrun{draw.default(filter.number=10, family="DaubLeAsymm")}
#
# Wow. What a great picture!
#
# Now how about a one-dimensional Daubechies extremal-phase scaling function
# with N=2
#
\dontrun{draw.default(filter.number=2, family="DaubExPhase")}
#
# Excellent! Now how about a two-dimensional Daubechies least-asymmetric
# N=6 wavelet
#
# N.b. we'll also reduce the resolution down a bit. If we used the default
# resolution of 8192 this would be probably too much for most computers.
#
\dontrun{draw.default(filter.number=6, family="DaubLeAsymm", dimension=2, res=256)}
#
# What a pretty picture!
}
\keyword{hplot}
\author{G P Nason}
