\name{putD.wd}
\alias{putD.wd}
\title{Puts a whole resolution level of mother wavelet coeffients into wd wavelet object.}
\description{
Makes a copy of the \code{\link{wd}} object, replaces some mother wavelet coefficients data in the copy, and then returns the copy. 
}
\usage{
\method{putD}{wd}(wd, level, v, boundary=FALSE, index=FALSE, \dots)
}
\arguments{
\item{wd}{Wavelet decomposition object into which you wish to insert the mother wavelet coefficients.}
\item{level}{the resolution level at which you wish to replace the mother wavelet coefficients.} 
\item{v}{the replacement data, this should be of the correct length.} 
\item{boundary}{If \code{boundary} is \code{FALSE} then only "real" data is replaced. If boundary is \code{TRUE} then the boundary (bookeeping) elements are replaced as well. Information about the lengths of the vectors can be found in the \code{\link{first.last}} database function and Nason and Silverman, 1994.}
\item{index}{If index is \code{TRUE} then the index numbers into the 1D array where the coefficient insertion would take place are returned. If index is \code{FALSE} (default) then the modified \code{wavelet decomposition} object is returned.}
\item{\dots}{any other arguments}
}
\details{
The function \code{\link{accessD}} obtains the mother wavelet coefficients for a particular level. The function \code{putD.wd} replaces father wavelet coefficients at a particular resolution level and returns a modified wd object reflecting the change.

The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as a linear vector. \code{PutD.wd} obtains information about where the smoothed data appears from the \code{fl.dbase} component of an \code{\link{wd.object}}, in particular the array 

\code{fl.dbase$first.last.d}

which gives a complete specification of index numbers and offsets for 

\code{wd.object$D}.

Note that this function is method for the generic function \code{\link{putD}}. When the \code{\link{wd.object}} is definitely a wd class object then you only need use the generic version of this function. 

Note also that this function only puts information into \code{\link{wd}} class objects. To extract coefficients from a \code{\link{wd}} object you have to use the \code{\link{accessD}} function (or more precisely, the \code{\link{accessD.wd}} method).
}
\value{
A \code{\link{wd}} class object containing the modified mother wavelet coefficients. 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{putD}}, \code{\link{wd.object}}, \code{\link{wd}}, \code{\link{accessD}},\code{\link{putD}}, \code{\link{first.last}}, 
}
\examples{
#
# Generate an EMPTY wd object:
#
zero <- rep(0, 16)
zerowd <- wd(zero)
#
# Put some random father wavelet coefficients into the object at
# resolution level 2. For the decimated wavelet transform there
# are always 2^i coefficients at resolution level i. So we have to
# insert 4 coefficients
#
mod.zerowd <- putD( zerowd, level=2, v=rnorm(4))
#
# If you plot mod.zerowd you will see that there are only 
# coefficients at resolution level 2 where you just put the coefficients.
#
# Now, for a time-ordered non-decimated wavelet transform object the
# procedure is exactly the same EXCEPT that there are going to be
# 16 coefficients at each resolution level. I.e.
#
# Create empty TIME-ORDERED NON-DECIMATED wavelet transform object
#
zerowdS <- wd(zero, type="station")
#
# Now insert 16 random coefficients at resolution level 2
#
mod.zerowdS <- putD(zerowdS, level=2, v=rnorm(16))
#
# Once more if you plot mod.zerowdS then there will only be
# coefficients at resolution level 2.
}
\keyword{manip}
\author{G P Nason}
