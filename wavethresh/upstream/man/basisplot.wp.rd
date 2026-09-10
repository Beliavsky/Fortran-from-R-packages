\name{basisplot.wp}
\alias{basisplot.wp}
\title{Function to graphically select a wavelet packet basis}
\usage{
\method{basisplot}{wp}(x, draw.mode=FALSE, \dots)
}
\arguments{
\item{x}{The \code{\link{wp.object}} for which you wish to select
	a basis graphically for.}
\item{draw.mode}{If TRUE then TWO graphics windows have to be open.
	Every time a packet is selected in the packet selection window,
	a representation of the wavelet packet basis function is drawn
	in the other window}
\item{\dots}{Other arguments}
}

\description{
Note, one or two (depending on the state of \code{draw.mode}) graphics
windows with mouse-clickable interfaces have to open to use this function.

Graphically select a wavelet packet basis associated with a 
wavelet packet object. Left-click selects packets, right click
exits the routine.
}

\details{
A wavelet packet basis described in WaveThresh using the node vector
object (class from \code{\link{MaNoVe.wp}}) which for wavelet packets
is \code{nvwp}. This function takes a \code{\link{wp.object}} object
and graphically depicts all possible basis function locations. The user
is then invited to click on different packets, these change colour.
When finished, the user right clicks on the graphic and the selected
basis is returned.

\emph{Note that the routine does not check to see whether the basis is
legal. You have to do this.} A legal basis can select packets from
different levels, however you can't select packets that both cover the
same packet index, however every packet index has to be covered.

A better function \emph{would} check basis legality!
}
\value{
An object of class \code{nvwp} which contains the specification
for the basis. 
}
\seealso{\code{\link{addpkt}}, \code{\link{InvBasis}}, \code{\link{MaNoVe.wp}}, \code{\link{plotpkt}},
	\code{\link{wp}}}
\author{G P Nason}
\keyword{hplot}
