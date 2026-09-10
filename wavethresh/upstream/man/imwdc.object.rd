\name{imwdc.object}
\alias{imwdc.object}
\title{Two-dimensional compressed wavelet decomposition objects.}
\description{
These are objects of classes 

\code{imwdc}

They represent a decomposition of an image with respect to a two-dimensional wavelet basis 
}
\details{
In previous releases the original image was stored as the "original" component of a imwd object. This is not done now as the resulting objects were excessively large. 

To uncompress this class of object back into an object of class \code{\link{imwd.object}} use the \code{\link{uncompress.imwdc}} function. 
}
\value{
The following components must be included in a legitimate `imwdc' object. 

\item{nlevelsWT}{number of levels in wavelet decomposition. If you raise 2 to the power of nlevels then you get the dimension of the image that you originally started with. }
\item{type}{If \code{type="wavelet"} then the image was decomposed according to the 2D Mallat pyramidal algorithm. If \code{type="station"} then the image was decomposed using the 2D spatially ordered non-decimated wavelet transform.}
\item{fl.dbase}{The first last database associated with the decomposition. For images, this list is not very useful as each level's components is stored as a list component, rather than being packaged up in a single vector as in the 1D case. Nevertheless the internals still need to know about fl.dbase to get the computations correct. See the help for \code{\link{first.last}} if you are a masochist. }
\item{filter}{A filter object as returned by the \code{\link{filter.select}} function. This component records the filter used in the decomposition. The reconstruction routines use this component to find out what filter to use in reconstruction. }
\item{wNLx}{The object will probably contain many components with names of this form. These are all the wavelet coefficients of the decomposition. In "wNLx" the "N" refers to the level number and the "x" refers to the direction of the coefficients with "1" being horizontal, "2" being vertical and "3" being diagonal. Note that imwdc objects do not contain scaling function coefficients. This would negate the point of having a compressed object.

Each vector stores its coefficients using an object of class compressed, i.e. the vector is run-length encoded on zeroes.

Note that the levels should be in numerically decreasing order, so if nlevelsWT is 5, then there will be w5L1, w5L2, w5L3 first, then down to w1L1, w1L2, and w1L3. Note that these coefficients store their data according to the \code{\link{first.last}} database \code{fl.dbase$first.last.d}, so refer to them using this. 

Note that if \code{type="wavelet"} then images at level N are subimages of side length \code{2^N} pixels. If the type component is \code{"station"} then each coefficient subimage is of the same dimension as the input image used to create this object.}
\item{w0Lconstant}{This is the coefficient of the bottom level scaling function coefficient. So for examples, if you used Haar wavelets this would be the sample mean of the data (scaled by some factor depending on the number of levels, nlevelsWT).}
\item{bc}{This component details how the boundaries were treated in the decomposition.} 
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{threshold.imwd}} function to represent a thresholded two-dimensional wavelet decomposition of a function. Some other functions return an object of class imwdc.  
}
\section{METHODS}{
The imwd class of objects has methods for the following generic functions: \code{\link{draw}}, \code{\link{imwr}}, \code{\link{nullevels}}, \code{\link{plot}}, \code{\link{print}}, \code{\link{summary}}, \code{\link{threshold.imwdc}}. 
}
\section{RELEASE}{
Version 3.5.3 Copyright Guy Nason 1994 
}
\seealso{
\code{\link{imwd}} \code{\link{imwd.object}}, \code{\link{threshold.imwd}}, \code{\link{uncompress.imwdc}}. 
}
\examples{
#
# Perform the standard two-dimensional DWT
# on the lennon image.
#
data(lennon)

lwd <- imwd(lennon)
#
# Now let's see how many horizontal detail coefficients there are at
# scale 6
#
length(lwd$w6L1)
# [1] 4096
#
# So the horizontal detail ``image'' at scale contains 64x64=4096 coefficients.
# A lot!
#
# Now, suppose we threshold this
# two-dimensional wavelet decomposition object
#
lwdT <- threshold(lwd)
#
# First of all. What is the class of the detail coefficients now? 
#
class(lwdT$w6L1)
# [1] "compressed"
#
# Aha. So this set of coefficients got compressed using the
# compress.default function.
#
# How many coefficients are being stored here?
#
lwdT$w6L1
# $position:
# [1]  173 2829 2832 2846
#
# $values:
# [1]  141.5455 -190.2810 -194.5714 -177.1791
# 
# $original.length:
# [1] 4096
#
# attr(, "class"):
# [1] "compressed"
#
# Wow! Only 4 coefficients are not zero. Wicked compression!
}
\keyword{classes}
\keyword{smooth}
\author{G P Nason}
