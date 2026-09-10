\name{wavegrow}
\alias{wavegrow}
\title{Interactive graphical tool to grow a wavelet synthesis}
\usage{
wavegrow(n = 64, filter.number = 10, family = "DaubLeAsymm", type = "wavelet",
    random = TRUE, read.value = TRUE, restart = FALSE)
}
\arguments{
\item{n}{Number of points in the decomposition}
\item{filter.number}{The wavelet filter.number to use,
	see \code{\link{filter.select}}}
\item{family}{The wavelet family to use in the reconstruction}
\item{type}{If \code{"wavelet"} then carry out the regular wavelet
	transform, otherwise if \code{"station"} do the nondecimated
	transform.}
\item{random}{If \code{TRUE} then iid Gaussian coefficients are inserted
	into the tableaux. If \code{FALSE} and \code{read.value=TRUE}
	then the user is promoted for a value, otherwise the value
	1 is inserted into the tableaux at the selected point.}
\item{read.value}{If \code{TRUE} then a value is read and used to insert
	that size of wavelet coefficient at the selected point. If \code{FALSE}
	then a coefficient of size 1 is inserted.}
\item{restart}{If \code{TRUE} then after a coefficient has been inserted,
	and plots done, the next selection causes all the coefficients
	to be reset to zero and a single coefficient inserted. This actually
	has the overall action of being able to select a coefficient
	location and view the size and shape of the wavelet produced.}
}
\description{
Use mouse to select which wavelets to enter a wavelet synthesis,
continually plot the reconstruction and the wavelet tableaux.
}
\details{
This function can perform many slightly different actions. However,
the basic idea is for a tableaux of wavelet coefficients to be displayed
in one graphics window, and the reconstruction of those coefficients
to be displayed in another graphics window.

Hence, two graphics windows, capable of plotting and mouse interaction
(e.g. X11, windows or quartz) with the locator function, are required to
be active.

When the function starts up an initial random tableaux is displayed and
its reconstruction.

The next step is for the user to select coefficients on the tableaux.
What happens next specifically depends on the arguments above. By default
selecting a coefficient causes that coefficient scale and location to
be identified, then a random sample is taken from a N(0,1) random variable
and assigned to that coefficient. Hence, the tableaux is updated, the
reconstruction with the new coefficient computed and both are plotted.

If \code{type="wavelet"} is used then decimated wavelets are used,
if \code{type="station"} then the time-ordered non-decimated wavelets
are used.

If \code{random=FALSE} then new values for the coefficients are either
selected (by asking the user for input) if \code{read.value=TRUE} or
the value of 1 is input. 

If \code{restart=TRUE} then the function merely displays the wavelet
associated with the selected coefficient. Hence, this option is useful
to demonstrate to people how wavelets from different points of the
tableaux have different sizes, scales and locations.

If the mouse locator function is exited (this can be a right-click in some
windowing systems, or pressing ESCAPE) then the function asks whether the
user wishes to continue. If not then the function returns the current
tableux. Hence, this function can be useful for users to build their
own tabeleaux.
}
\value{
The final tableaux.
}
\seealso{\code{\link{wd}}}
\author{G P Nason}
\keyword{hplot}
\keyword{iplot}
