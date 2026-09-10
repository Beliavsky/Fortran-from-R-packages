\name{accessD.wd}
\alias{accessD.wd}
\title{Get detail (mother wavelet) coefficients from wavelet object (wd).}
\description{
This function extracts and returns a vector of mother wavelet coefficients, corresponding to a particular resolution level, from a \code{\link{wd}} wavelet decomposition object.
 
The pyramid of coefficients in a wavelet decomposition (returned from the \code{\link{wd}} function, say) are packed into a single vector in WaveThresh. 
}
\usage{
\method{accessD}{wd}(wd, level, boundary=FALSE, aspect="Identity", \dots)
}
\arguments{
\item{wd}{Wavelet decomposition object from which you wish to extract the mother wavelet coefficients.}
\item{level}{The resolution level at which you wish to extract coefficients.}
\item{boundary}{some methods of wavelet transform computation handle the boundaries by keeping some extra bookkeeping coefficients at either end of a resolution level. If this argument is TRUE then these bookkeeping coefficients are returned when the mother wavelets are returned. Otherwise, if FALSE, these coefficients are not returned.}
\item{aspect}{The aspect argument permits the user to supply a function to modify the returned coefficients. The function is applied to the vector of coefficients before it is returned. This can be useful, say, with the complex DWT where you could supply aspect="Mod" if you wanted to return the modulus of the coefficients at a given resolution level. The default argument, "Identity", ensures that the coefficients are not modified before returning.}
\item{\dots}{any other arguments}
}
\details{
The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as a linear vector. AccessD obtains information about where the smoothed data appears from the \code{fl.dbase} component of an \code{\link{wd}} object, in particular the array
 
\code{fl.dbase$first.last.d}

which gives a complete specification of index numbers and offsets for
 
\code{wd.object$D}.

Note that this function is a method for the generic function \code{\link{accessD}}. 

Note also that this function only retrieves information from \code{\link{wd}} class objects. To insert coefficients into \code{\link{wd}} objects you have to use the \code{\link{putD}} function (or more precisely, the \code{\link{putD.wd}} method).
}
\value{
A vector containing the mother wavelet coefficients at the required resolution level (the coefficients might have been modified depending on the value of the aspect argument).
}

\references{
Mallat, S. G. (1989) A theory for multiresolution signal decomposition: the wavelet representation. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}. \bold{11}, 674--693. 

Nason, G. P. and Silverman, B. W. (1994). The discrete wavelet transform in S. \emph{Journal of Computational and Graphical Statistics}, \bold{3}, 163--191
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wr}}, \code{\link{wd}}, \code{\link{accessD}}, \code{\link{filter.select}}, \code{\link{threshold}}
}
\examples{
#
# Get the 4th resolution level of wavelet coefficients. 
#
dat <- rnorm(128)
accessD(wd(dat), level=4)
}
\keyword{manip}
\author{G P Nason}
