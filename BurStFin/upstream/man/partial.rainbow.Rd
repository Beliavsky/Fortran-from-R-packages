\name{partial.rainbow}
\alias{partial.rainbow}
\title{ Create Palette Function for Part of Rainbow}
\description{
Returns a function suitable as the \code{color.palette} argument
to \code{filled.contour} that contains a specified portion of the
rainbow.
}
\usage{
partial.rainbow(start = 0, end = 0.35)
}
\arguments{
  \item{start}{ 
	a number giving where the colors should start.
	Valid numbers range from 0 (red) to 1 (also red).
	}
  \item{end}{ 
	a number giving where the colors should end.
	Valid numbers range from 0 (red) to 1 (also red).
	}
}
\value{
a function similar to \code{rainbow} but with the \code{start}
and \code{end} arguments (possibly) changed.
}
\section{Details }{
This function was made to facilitate the construction of contour
plots in the \code{alpha.proxy} function, but is of general use.
}
\references{ 
Burns, Patrick (2003). "Portfolio Sharpening". Working Paper,
Burns Statistics http://www.burns-stat.com/.
	}
\section{Revision }{
This help was last revised 2010 January 05.
}

\seealso{ 
\code{\link{rainbow}}, \code{\link{filled.contour}}, \code{\link{alpha.proxy}}.
}
\examples{
ap1 <- alpha.proxy(cor.man=seq(0, .3, len=50), 
	vol.man=seq(.15, .25, leng=60), plot=FALSE)
filled.contour(ap1, color.palette=partial.rainbow(start=.05, end=.3))
}
\keyword{color}
