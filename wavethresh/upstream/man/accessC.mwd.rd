\name{accessC.mwd}
\alias{accessC.mwd}
\title{Get Smoothed Data from Wavelet Structure}
\description{
  The smoothed and original data from a multiple wavelet decomposition structure, \code{\link{mwd.object}}  ect  (e.g. returned from \code{\link{mwd}}) are packed into a single matrix in that structure.  TRUE  his function extracts the data corresponding to a particular resolution level.
}
\usage{
\method{accessC}{mwd}(mwd, level = nlevelsWT(mwd), \dots)
}
\arguments{
   \item{mwd}{Multiple wavelet decomposition structure from which you wish to extract the smoothed or original data if the structure is from a wavelet decomposition, or the reconstructed data if the structure is from a wavelet reconstruction.}
   \item{level}{The level that you wish to extract. By default, this is the level with most detail (in the case of structures from a decomposition this is the original data, in the case of structures from a reconstruction this is the top-level reconstruction).}
\item{\dots}{any other arguments}
}
\details{
The \link{mwd} function produces a wavelet decomposition structure.

For decomposition, the top level contains the original data, and subsequent lower levels contain the successively smoothed data. So if there are \code{mwd$filter$npsi*2^m} original data points (\code{mwd$filter$npsi} is the multiplicity of wavelets), there will be \code{m+1} levels indexed 0,1,...,m. So 

\code{accessC.mwd(Mwd, level=m)}

pulls out the original data, as does

\code{accessC.mwd(mwd)}

To get hold of lower levels just specify the level that you're interested in, e.g. 

\code{accessC.mwd(mwd, level=2)}

Gets hold of the second level.

The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as a linear matrix of coefficients. AccessC obtains information about where the smoothed data appears from the fl.dbase component of mwd, in particular the array \code{fl.dbase$first.last.c} which gives a complete specification of index numbers and offsets for \code{mwd$C}. 

Note also that this function only gets information from \link{mwd} class objects. To \emph{put} coefficients into \link{mwd} structures you have to use the \link{putC.mwd} function. 

See Downie and Silverman, 1998. 
}
\value{
A matrix with \code{mwd$filter$npsi} rows containing the extracted data of all the coefficients at that level.
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6.)}
\seealso{ 
\code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}
 }
\examples{
#
# Get the 3rd level of smoothed data from a decomposition
#
dat <- rnorm(32)
accessC.mwd(mwd(dat), level=3)
}
\keyword{manip}
\author{G P Nason}
