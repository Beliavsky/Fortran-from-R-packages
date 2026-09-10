\name{imwd.object}
\alias{imwd.object}
\title{Two-dimensional wavelet decomposition objects.}
\description{
These are objects of classes 

\code{imwd}

They represent a decomposition of an image with respect to a two-dimensional wavelet basis (or tight frame in the case of the two-dimensional (space-ordered) non-decimated wavelet decomposition). 
}
\details{
In previous releases the original image was stored as the "original" component of a imwd object. This is not done now as the resulting objects were excessively large. 
}
\value{
The following components must be included in a legitimate `imwd' object. 

\item{nlevelsWT}{number of levels in wavelet decomposition. If you raise 2 to the power of nlevels then you get the dimension of the image that you originally started with. }
\item{type}{If \code{type="wavelet"} then the image was decomposed according to the 2D Mallat pyramidal algorithm. If \code{type="station"} then the image was decomposed using the 2D spatially ordered non-decimated wavelet transform.}
\item{fl.dbase}{The first last database associated with the decomposition. For images, this list is not very useful as each level's components is stored as a list component, rather than being packaged up in a single vector as in the 1D case. Nevertheless the internals still need to know about fl.dbase to get the computations correct. See the help for \code{\link{first.last}} if you are a masochist. }
\item{filter}{A filter object as returned by the \code{\link{filter.select}} function. This component records the filter used in the decomposition. The reconstruction routines use this component to find out what filter to use in reconstruction. }
\item{wNLx}{The object will probably contain many components with names of this form. These are all the wavelet coefficients of the decomposition. In "wNLx" the "N" refers to the level number and the "x" refers to the direction of the coefficients with "1" being horizontal, "2" being vertical and "3" being diagonal and "4" corresonding to scaling function coefficients at the given resolution level. Note that the levels should be in numerically decreasing order, so if nlevelsWT is 5, then there will be w5L1, w5L2, w5L3 first, then down to w1L1, w1L2, and w1L3. Note that these coefficients store their data according to the \code{\link{first.last}} database \code{fl.dbase$first.last.d}, so refer to them using this. 
Note that if \code{type="wavelet"} then images at level N are subimages of side length \code{2^N} pixels. If the type component is \code{"station"} then each coefficient subimage is of the same dimension as the input image used to create this object.}
\item{w0Lconstant}{This is the coefficient of the bottom level scaling function coefficient. So for examples, if you used Haar wavelets this would be the sample mean of the data (scaled by some factor depending on the number of levels, nlevelsWT).}
\item{bc}{This component details how the boundaries were treated in the decomposition.} 
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{imwd}} function to represent a two-dimensional (possibly space-ordered non-decimated) wavelet decomposition of a function. Many other functions return an object of class imwd. 
}
\section{METHODS}{
The imwd class of objects has methods for the following generic functions: \code{\link{compress}}, \code{\link{draw}}, \code{\link{imwr}}, \code{\link{nullevels.imwd}}, \code{\link{plot}}, \code{\link{print}}, \code{\link{summary}}, \code{\link{threshold.imwd}}. 
}
\section{RELEASE}{
Version 3.5.3 Copyright Guy Nason 1994 
}
\seealso{
\code{\link{imwd}}
}
\keyword{classes}
\keyword{smooth}
\author{G P Nason}
