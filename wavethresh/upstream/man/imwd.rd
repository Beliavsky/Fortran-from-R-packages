\name{imwd}
\alias{imwd}
\title{Two-dimensional wavelet transform (decomposition). }
\description{
This function can perform two types of two-dimensional discrete wavelet transform (DWT). The standard transform (\code{type="wavelet"}) computes the 2D DWT according to Mallat's pyramidal algorithm (Mallat, 1989). The spatially ordered non-decimated 2D DWT (NDWT) (\code{type="station"}) contains all possible spatially shifted versions of the DWT. The order of computation of the DWT is O(n), and it is O(n log n) for the NDWT if n is the number of pixels.
}
\usage{
imwd(image, filter.number=10, family="DaubLeAsymm", type="wavelet",
bc="periodic", RetFather=TRUE, verbose=FALSE)

}
\arguments{
\item{image}{A square matrix containing the image data you wish to decompose. The sidelength of this matrix must be a power of 2.} 
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. The options are "DaubExPhase" and "DaubLeAsymm".}
\item{type}{specifies the type of wavelet transform. This can be "wavelet" (default) in which case the standard 2D DWT is performed (as in previous releases of WaveThresh). If type is "station" then the 2D spatially-ordered non-decimated DWT is performed. At present, only periodic boundary conditions can be used with the 2D spatially ordered non-decimated wavelet transform.}
\item{bc}{specifies the boundary handling. If bc=="periodic" the default, then the function you decompose is assumed to be periodic on it's interval of definition, if bc=="symmetric" then the function beyond its boundaries is assumed to be a symmetric reflection of the function in the boundary. The symmetric option was the implicit default in releases prior to 2.2. Note that only periodic boundary conditions are valid for the 2D spatially-ordered non-decimated wavelet transform.}
\item{RetFather}{If \code{TRUE} then this argument causes the scaling function coefficients at each resolution level to be returned as well as the wavelet coefficients. If \code{FALSE} then no scaling function coefficients are returned. The opportunity of returning father wavelet coefficients has been added since previous versions of WaveThresh.}
\item{verbose}{Controls the printing of "informative" messages whilst the computations progress. Such messages are generally annoying so it is turned off by default.}
}
\details{
The 2D algorithm is essentially the application of many 1D filters. First, the columns are attacked with the smoothing (H) and bandpass (G) filters, and the rows of each of these resultant images are attacked again with each of G and H, this results in 4 images. Three of them, GG, GH, and HG correspond to the highest resolution wavelet coefficients. The HH image is a smoothed version of the original and can be further attacked in exactly the same way as the original image to obtain GG(HH), GH(HH), and HG(HH), the wavelet coefficients at the second highest resolution level and HH(HH) the twice-smoothed image, which then goes on to be further attacked. 

If \code{RetFather=TRUE} then the results of the HH smooth (the scaling function coefficients) are returned additionally. 

There are now two methods of handling "boundary problems". If you know that your function is periodic (on it's interval) then use the bc="periodic" option, if you think that the function is symmetric reflection about each boundary then use bc="symmetric". If you don't know then it is wise to experiment with both methods, in any case, if you don't have very much data don't infer too much about your decomposition! If you have loads of data then don't worry too much about the boundaries. It can be easier to interpret the wavelet coefficients from a bc="periodic" decomposition, so that is now the default. 

The spatially-ordered non-decimated DWT contains all spatial (toroidal circular) shifts of the standard DWT. 

The standard DWT is orthogonal, the spatially-ordered non-decimated transform is most definitely not. This has the added disadvantage that non-decimated wavelet coefficients, even if you supply independent normal noise. This is unlike the standard DWT where the coefficients are independent (normal noise). 

The two-dimensional packet-ordered non-decimated discrete wavelet transform is computed by the \code{\link{wst2D}} function. 

}
\value{
An object of class \code{\link{imwd.object}} containing the two-dimensional wavelet transform (possibly spatially-ordered non-decimated). 
}
\section{RELEASE}{Version 3.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wd}}, \code{\link{imwd.object}}, \code{\link{filter.select}} 
}
\examples{
data(lennon)
#
# Let's use the lennon test image
#
\dontrun{image(lennon)}
#
# Now let's do the 2D discrete wavelet transform
#
lwd <- imwd(lennon)
#
# Let's look at the coefficients
#
\dontrun{plot(lwd)}
}
\author{G P Nason}
\keyword{smooth}
