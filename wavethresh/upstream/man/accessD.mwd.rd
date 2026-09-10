\name{accessD.mwd}
\alias{accessD.mwd}
\title{Get wavelet coefficients from multiple wavelet structure (mwd).}
\description{
 The wavelet coefficients from a multiple wavelet decomposition structure, \code{\link{mwd.object}}, (e.g. returned from \code{\link{mwd}}) are packed into a single matrix in that structure. This function extracts the coefficients corresponding to a particular resolution level.
}
\usage{
\method{accessD}{mwd}(mwd, level, \dots)
}
\arguments{
 \item{mwd}{Multiple wavelet decomposition structure from which you wish to extract the expansion coefficients.}
\item{level}{The level that you wish to extract. If the "original" data has \code{mwd$filter$npsi*2^m} data points (\code{mwd$filter$npsi} being the multiplicity of the multiple wavelets) then there are m possible levels that you could want to access, indexed by 0,1,...,(m-1) }
\item{\dots}{any other arguments}
}
\details{
The \code{\link{mwd}} function produces a \code{multiple wavelet decomposition object} .
 
The need for this function is a consequence of the pyramidal structure of
\code{Mallats algorithm} and the memory efficiency gain achieved by storing
the pyramid as a linear matrix.
AccessD obtains information about where the coefficients appear from the
fl.dbase component of \code{\link{mwd}},
in particular the array \code{fl.dbase$first.last.d} which gives a complete
specification of index numbers and offsets for \code{mwd$D}. 

Note that this function and \code{\link{accessC}} only work on objects of class \code{\link{mwd}} to \emph{extract} coefficients. You have to use
\code{\link{putD.mwd}} to insert wavelet coefficients into a \code{\link{mwd}} object. 

See Downie and Silverman, 1998. 
}
\value{
A matrix with \code{mwd$filter$npsi} rows containing the extracted coefficients.
}

\section{RELEASE}{Tim Downie 1995-6}
\seealso{
\code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}
}
\examples{
#
# Get the 3rd level of smoothed data from a decomposition
#
data(ipd)
accessD.mwd(mwd(ipd), level=3)
}
\keyword{manip}
\author{G P Nason}
