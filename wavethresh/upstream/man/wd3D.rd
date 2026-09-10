\name{wd3D}
\alias{wd3D}
\title{Three-dimensional discrete wavelet transform}
\description{
This function performs the 3D version of Mallat's discrete wavelet transform (see Mallat, 1989, although this paper does not describe in detail the 3D version the extension is trivial). The function assumes
\emph{periodic} boundary conditions. 
}
\usage{
wd3D(a, filter.number=10, family="DaubLeAsymm")
}
\arguments{
\item{a}{A three-dimensional array constructed using the S-Plus \code{array()} function. Each dimension of the array should be equal to the same power of two.}
\item{filter.number}{This selects the smoothness of wavelet that you want to use in the decomposition. By default this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments.}
\item{family}{specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for \code{\link{filter.select}} for more possibilities.} 
}
\details{
This function implements a straightforward extension of Mallat's, (1989) one- and two-dimensional DWT. The algorithm recursively applies all possible combinations of the G and H detail and smoothing filters to each of the dimensions thus forming 8 different sub-blocks which we label HHH, GHH, HGH, GGH, HHG, GHG, HGG, and GGG. The algorithm recurses on the HHH component of each level (these are the father wavelet coefficients). 

Making an analogy to the 2D transform where HH, HG, HG and GG is produced at each resolution level: the HG and GH correspond to "horizontal" and "vertical" detail and GG corresponds to "diagonal detail". The GGG corresponds to the 3D "diagonal" version, HGG corresponds to smoothing in dimension 1 and "diagonal" detail in dimensions 2 and 3, and so on. I don't think there are words in the English language which adequately describe "diagonal" in 3D --- maybe cross detail? 
}
\value{
An object of class \code{\link{wd3D}}. 
}
\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1997}
\seealso{
\code{\link{wd}}, \code{\link{imwd}}, \code{\link{accessD.wd3D}}, \code{\link{print.wd3D}}, \code{\link{putD.wd3D}}, \code{\link{putDwd3Dcheck}}, \code{\link{summary.wd3D}}, \code{\link{threshold.wd3D}}, \code{\link{wd3D.object}}, \code{\link{wr3D}}. 
}
\examples{
#
# Generate some test data: 512 standard normal observations in an 8x8x8
# array.
#
test.data.3D <- array(rnorm(8*8*8), dim=c(8,8,8))
#
# Now do the 3D wavelet transform
#
tdwd3D <- wd3D(test.data.3D)
#
# See examples explaining  the 3D wavelet transform.
#
}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
