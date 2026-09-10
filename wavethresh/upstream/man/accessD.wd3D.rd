\name{accessD.wd3D}
\alias{accessD.wd3D}
\title{Get wavelet coefficients from 3D wavelet object}
\description{
This function extracts and returns arrays of wavelet coefficients, corresponding to a particular resolution level, from a \code{\link{wd}} wavelet decomposition object.
 
The pyramid of coefficients in a wavelet decomposition (returned from the \code{\link{wd3D}} function, say) are packed into a single array in \code{WaveThresh3}.
}
\usage{
\method{accessD}{wd3D}(obj, level = nlevelsWT(obj)-1, block, \dots)
} 
\arguments{
\item{obj}{3D Wavelet decomposition object from which you wish to extract the wavelet coefficients.}
\item{level}{The resolution level at which you wish to extract coefficients. The minimum level you can enter is 0, the largest is one less than the number of nlevelsWT stored in the obj object.}
\item{block}{if block is missing then a list containing all of the wavelet coefficient blocks GGG, GGH, GHG, GHH, HGG, HGH, HHG (and HHH, if level=0) is returned. Otherwise block should be one of the character strings GGG, GGH, GHG, GHH, HGG, HGH, HHG and then only that sub-block is returned from the resolution level specified.}
\item{\dots}{any other arguments}
}
\details{
The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as a array. 

Note that this functiOn is a method for the generic function \code{\link{accessD}}.
}
\value{
If the block is missing then a list is returned containing all the sub-blocks of coefficients for the specificed resolution \code{level}.
 
Otherwise the block character string specifies which sub-block of coefficients to return. 
}

\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1997 }
\seealso{
\code{link{accessD}}, \code{link{print.wd3D}}, \code{link{putD.wd3D}},  \code{link{putDwd3Dcheck}},  \code{link{summary.wd3D}},  \code{link{threshold.wd3D}},  \code{link{wd3D}},  \code{link{wd3D object}},  \code{link{wr3D}}.
}
\examples{
#
# Generate some test data
#
a <- array(rnorm(8*8*8), dim=c(8,8,8))
#
# Perform the 3D DWT
#
awd3D <- wd3D(a)
#
# How many levels does this object have?
#
nlevelsWT(awd3D)
# [1] 3
#
# So conceivably we could access levels 0, 1 or 2.
#
# Ok. Let's get the level 1 HGH sub-block coefficients:
#
accessD(awd3D, level=1, block="HGH")
#
#, , 1
#           [,1]      [,2] 
#[1,]  0.8359289 1.3596832
#[2,] -0.1771688 0.2987303
#
#, , 2
#           [,1]       [,2] 
#[1,] -1.2633313 1.00221652
#[2,] -0.3004413 0.04728019
#
# This was a 3D array of dimension size 2 (8 -> 4 -> 2, level 3, 2 and then 1)
#
#
# Let's do the same call except this time don't specify the block arg.
#
alllev1 <- accessD(awd3D, level=1)
#
# This new object should be a list containing all the subblocks at this level.
# What are the components?
#
names(alllev1)
#[1] "GHH" "HGH" "GGH" "HHG" "GHG" "HGG" "GGG"
#
# O.k. Let's look at HGH again
#
alllev1$HGH
#
#, , 1
#           [,1]      [,2] 
#[1,]  0.8359289 1.3596832
#[2,] -0.1771688 0.2987303
#
#, , 2
#           [,1]       [,2] 
#[1,] -1.2633313 1.00221652
#[2,] -0.3004413 0.04728019
#
# Same as before.
#
}
\keyword{manip}
\author{G P Nason}
