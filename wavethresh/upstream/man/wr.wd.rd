\name{wr.wd}
\alias{wr.wd}
\title{Wavelet reconstruction for wd class objects (inverse discrete wavelet transform). }
\description{
This function performs the reconstruction stage of Mallat's pyramid algorithm (Mallat 1989), i.e. the discrete inverse wavelet transform. The actual transform is performed by some C code, this is dynamically linked into S (if your machine can do this). 
}
\usage{
\method{wr}{wd}(wd, start.level = 0, verbose = FALSE, bc = wd$bc,
	return.object = FALSE, filter.number = wd$filter$filter.number,
	family = wd$filter$family, \dots)
}
\arguments{
\item{wd}{A wavelet decomposition object as returned by \code{\link{wd}}, and described in the help for that function and the help for \code{\link{wd.object}}.}
\item{start.level}{The level you wish to start reconstruction at. The is usually the first (level 0). This argument is ignored for a wd object computed using the ``wavelets on the interval'' transform (i.e. using the \code{bc="interval"} option of \code{\link{wd}}.}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
\item{bc}{The boundary conditions used. Usually these are determined by those used to create the supplied wd object, but you sometimes change them with possibly silly results.}
\item{filter.number}{The filter number of the wavelet used to do the reconstruction. Again, as for bc, you should probably leave this argument alone. Ignored if the bvc component of the \code{\link{wd}} object is "\code{interval}".} 
\item{family}{The type of wavelet used to do the reconstruction. You can change this argument from the default but it is probably NOT wise. Ignored if the bvc component of the \code{\link{wd}} object is "\code{interval}".} 
\item{return.object}{If this is F then the top level of the reconstruction is returned (this is the reconstructed function at the highest resolution). Otherwise if it is T the whole wd reconstructed object is returned. Ignored if the \code{bvc} component of the \code{\link{wd}} object is "\code{interval}".} 
\item{\dots}{any other arguments}
}
\details{
The code implements Mallat's inverse pyramid algorithm. In the reconstruction the quadrature mirror filters G and H are supplied with c0 and d0, d1, ... d(m- 1) (the wavelet coefficients) and rebuild c1,..., cm. 

If the \code{bc} component of the \code{\link{wd}} object is "\code{interval}" then the \code{wr.int} function which implements the inverse ``wavelet on the interval'' transform due to Cohen, Daubechies and Vial, 1993 is used instead. 
}
\value{
Either a vector containing the top level reconstruction or an object of class wd containing the results of the reconstruction, details to be found in help for \code{\link{wd.object}}. 
}
\section{RELEASE}{Version 3 Copyright Guy Nason 1994
Integration of ``wavelets on the interval'' code by Piotr Fryzlewicz and Markus Monnerjahn was at Version 3.9.6, 1999. }
\seealso{
\code{\link{wd}}, \code{\link{wr.int}}, \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{filter.select}}, \code{\link{plot.wd}}, \code{\link{threshold}} 
}
\examples{
#
# Take the wd object generated in the examples to wd (called wds)
#
# Invert this wd object
#
#yans <- wr(wds)
#
# Compare it to the original, called y
#
#sum((yans-y)^2)
#[1] 9.805676e-017
#
# A very small number
#
}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
