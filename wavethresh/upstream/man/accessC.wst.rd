\name{accessC.wst}
\alias{accessC.wst}
\title{Get smoothed data from packet ordered non-decimated wavelet object (wst)}
\description{
  The smoothed data from a packet ordered non-decimated wavelet object (returned from \code{\link{wst}}) are stored in a matrix. This function extracts all the coefficients corresponding to a particular resolution level.
}
\usage{
\method{accessC}{wst}(wst, level, aspect, \dots)
}
\arguments{
\item{wst}{Packet ordered non-decimated wavelet object from which you wish to extract the smoothed or original data (if the object is directly from a packet ordered non-decimated wavelet transform of some data).}
\item{level}{The level that you wish to extract. This can range from zero (the coarsest coefficients) to nlevelsWT(wstobj) which returns the original data.}
\item{aspect}{Applies function to coefficients before return. Supplied
	as a character string which gets converted to a function.
	For example "Mod" which returns the absolute values of the coefficients}
\item{\dots}{Other arguments}
}
\value{
A vector of the extracted data. 
}
\details{
The \code{\link{wst}} function performs a packet-ordered non-decimated wavelet transform. This function extracts all the father wavelet coefficients at a particular resolution level specified by \code{level}.

Note that coefficients returned by this function are in \emph{packet order}. They can be used \emph{as is} but for many applications it might be more useful to deal with the coefficients in packets: see the function \code{\link{getpacket.wst}} for further details. 
}
\references{
Nason, G. P. and Silverman, B. W. (1994). The discrete wavelet transform in S. \emph{Journal of Computational and Graphical Statistics}, \bold{3}, 163--191. 
}
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{accessC}}, \code{\link{getpacket.wst}}
}
\examples{
#
# Get the 3rd level of smoothed data from a decomposition
#
dat <- rnorm(64)
accessC(wst(dat), level=3)
}
\keyword{manip}
\author{G P Nason}
