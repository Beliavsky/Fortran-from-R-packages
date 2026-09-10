\name{accessD.wst}
\alias{accessD.wst}
\title{Get mother wavelet coefficients from a packet ordered non-decimated wavelet object (wst).}
\description{
The mother wavelet coefficients from a packet ordered non-decimated wavelet object, \code{\link{wst}}, are stored in a matrix. This function extracts all the coefficients corresponding to a particular resolution level.
}
\usage{
\method{accessD}{wst}(wst, level, aspect = "Identity", \dots)
}
\arguments{
\item{wst}{Packet ordered non-decimated wavelet object from which you wish to extract the mother wavelet coefficients.}
\item{level}{The level that you wish to extract. This can range from zero (the coarsest coefficients) to nlevelsWT(wstobj) which returns the original data. }
\item{aspect}{Function to apply to coefficient before return.
	Supplied as a character argument which gets converted to a function.
	For example, "Mod" which returns the absolute value of complex-valued
	coefficients.}
\item{\dots}{Other arguments}
}
\details{
The \code{\link{wst}} function performs a packet-ordered non-decimated wavelet transform. This function extracts all the mother wavelet coefficients at a particular resolution level specified by \code{level}.
 
Note that coefficients returned by this function are in \emph{packet order}. They can be used \emph{as is} but for many applications it might be more useful to deal with the coefficients in packets: see the function \code{\link{getpacket.wst}} for further details. 

Note that all the coefficients here are those of mother wavelets. The non-decimated transform efficiently computes all possible shifts of the discrete wavelet transform computed by \code{\link{wd}}.
}
\value{
A vector of the extracted coefficients. 
} 
\references{
Nason, G.P. and Silverman, B.W. The stationary wavelet transform and some statistical applications. 
}
\seealso{
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{accessD}}, \code{\link{getpacket.wst}}
} 
\examples{
#
# Get the 4th level of mother wavelet coefficients from a decomposition
#
dat <- rnorm(128)
accessD(wst(dat), level=4)
}
\keyword{manip}
\author{G P Nason}
