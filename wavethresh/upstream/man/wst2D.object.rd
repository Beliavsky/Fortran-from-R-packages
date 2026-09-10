\name{wst2D.object}
\alias{wst2D.object}
\title{(Packet ordered) Two-dimensional nondecimated wavelet transform decomposition objects.}
\description{
These are objects of class \code{wst2D} They represent a decomposition of a function with respect to a set of (all possible) shifted two-dimensional wavelets. They are a 2D extension of the \code{\link{wst.object}}. 
}
\details{
To retain your sanity we recommend that the coefficients from a \code{wst2D} object be extracted or replaced using 


\itemize{
\item{\code{\link{getpacket.wst2D}} to obtain individual packets of either father or mother wavelet coefficients.}
\item{\code{\link{putpacket.wst2D}} to insert coefficients.}
}
 
You can obtain the coefficients directly from the \code{wst2D$wst2D} component but you have to understand their organization described above.}
\value{
The following components must be included in a legitimate \code{wst2D} object. 
\item{wst2D}{This a three-dimensional array. Suppose that the original image that created the \code{wst2D} object is n x n. Then the dimension of the \code{wst2D} array is [ nlevelsWT, 2n, 2n]. The first index of the array refers to the resolution level of the array with "resolution level = index - 1" (so, e.g. the coarsest scale detailed is stored at index 1 and the finest at index nlevels). For a given resolution level (selected first index) the associated 2n x 2n matrix contains the two-dimensional non-decimated wavelet coefficients for that level packed as follows. 

At the finest resolution level the 2n x 2n coefficient image may be broken up into four n x n subimages. Each of the four images corresponds to data shifts in the horizontal and vertical directions. The top left image corresponds to ``no shift'' and indeed the top left image corresponds to the coefficients obtained using the decimated 2D wavelet transform as obtained using the \code{\link{imwd}} function. The top right image corresponds to a horizontal data shift; the bottom left to a vertical data shift and the bottom right corresponds to both horizontal and vertical data shift. 

Within each of the four n x n images named in the previous paragraph are again 4 subimages each of dimension n/2 x n/2. These correspond to (starting at the top left and moving clockwise) the smooth (CC), horizontal detail (DC), diagonal detail (DD) and vertical detail (CD). 

At coarser resolution levels the coefficients are smaller submatrices corresponding to various levels of data shifts and types of detail (smooth, horizontal, vertical, diagonal). 

We strongly recommend the use of the \code{\link{getpacket.wst2D}} and \code{\link{putpacket.wst2D}} functions to remove and replace coefficients from \code{wst2D}} objects.
\item{nlevelsWT}{The number of levels in the decomposition. If you raise 2 to the power of 2 \code{nlevels} you get the number of data points used in the decomposition.}
\item{filter}{a list containing the details of the filter that did the decomposition (equivalent to the return value from the \code{\link{filter.select}} function).}
\item{date}{The date that the transform was performed or the \code{wst2D} was modified.} 
}
\section{GENERATION}{
This class of objects is returned from the \code{\link{wst2D}} function which computes the \emph{packets-ordered} two-dimensional non-decimated wavelet transform (effectively all possible shifts of the standard two-dimensional discrete wavelet transform). 

Many other functions return an object of class \code{wst2D}. 
}
\section{METHODS}{
The wst2D class of objects has methods for the following generic functions: \code{\link{AvBasis}}, \code{\link{getpacket}}. \code{\link{plot}}, \code{\link{print}}, \code{\link{putpacket}}, \code{\link{summary}}, 
}
\section{RELEASE}{
Version 3.5.3 Copyright Guy Nason 1994
}
\seealso{
\code{\link{wst2D}}
}
\keyword{classes}
\keyword{smooth}
\author{G P Nason}
