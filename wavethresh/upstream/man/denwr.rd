\name{denwr}
\alias{denwr}
\title{Wavelet reconstruction for density estimation.  }
\usage{
denwr(wd, start.level=0, verbose=FALSE, bc=wd$bc,
    return.object=FALSE, filter.number=wd$filter$filter.number,
    family=wd$filter$family)
}
\arguments{
\item{wd}{Wavelet decomposition object to reconstruct}
\item{start.level}{The level you wish to start the reconstruction at.
	This is usually the first level (level 0).
	Note that this option assumes the coarsest level is labelled 0,
	so it is best to think of this argument as
	"the number of levels up from the coarsest level to start the
	reconstruction".}
\item{verbose}{Controls the printing of "informative" messages whilst the
	computations progress. Such messages are generally annoying
	so it is turned off by default.}
\item{bc}{The boundary conditions used. These should be determined by those
	used to create the supplied \code{\link{wd.object}} object.
	In the case of density estimation they are "zero".}
\item{filter.number}{The filter number of the wavelet used to do the
	reconstruction. Again, as for bc, you should probably leave this argument alone.}
\item{family}{The type of wavelet used to do the reconstruction. You can change this argument from the default but it is probably NOT wise.}
\item{return.object}{If this is FALSE then the top level of the reconstruction
	is returned (this is the reconstructed function at the highest
	resolution). Otherwise, if it is TRUE, the whole wd
	reconstructed object is returned.}
}
\description{
Performs wavelet reconstruction for density estimation.
}
\details{
This is the same as \code{\link{wr.wd}},
except that it can handle zero boundary conditions.
}
\value{
Either a vector containing the top level reconstruction or an object of
	class \code{\link{wd.object}} containing the results of the reconstruction.
}
\author{David Herrick}
\keyword{smooth}
