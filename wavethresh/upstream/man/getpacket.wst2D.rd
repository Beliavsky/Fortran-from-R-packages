\name{getpacket.wst2D}
\alias{getpacket.wst2D}
\title{Get packet of coefficients from a two-dimensional non-decimated wavelet object (wst2D).}
\description{
This function extracts and returns a packet of coefficients from a two-dimensional non-decimated wavelet (\code{\link{wst2D}}) object. 
}
\usage{
\method{getpacket}{wst2D}(wst2D, level, index, type="S", Ccode=TRUE, \dots)
}
\arguments{
\item{wst2D}{2D non-decimated wavelet object from which you wish to extract a packet from.} 
\item{level}{The resolution level of the coefficients that you wish to extract. Can range from 0 to \code{\link{nlevelsWT}}(wpst)-1.} 
\item{index}{The index number within the resolution level of the packet of coefficients that you wish to extract. Index is a base-4 number which is r digits long. Each digit can be 0, 1, 2 or 3 corresponding to no shifts, horizontal shift, vertical shift or horizontal and vertical shifts. The number r indicates the depth of the resolution level from the data resolution i.e. where \code{r = nlevelsWT - level}.
 
Where there is a string of more than one digit the left most digits correspond to finest scale shift selection, the right most digits to the coarser scales (I think).}
\item{type}{This is a one letter character string: one of "S", "H", "V" or "D" for the smooth coefficients, horizontal, vertical or diagonal detail.} 
\item{Ccode}{If \code{T} then fast C code is used to obtain the packet, otherwise slow SPlus code is used. Unless you have some special reason always use the C code (and leave the argument at its default).} 
\item{\dots}{any other arguments}
}
\details{
The \code{\link{wst2D}} function creates a \code{\link{wst2D}} class object. Starting with a smooth the operators H, G, GS and HS (where G, H are the usual Mallat operators and S is the shift-by-one operator) are operated first on the rows and then the columns: i.e. so each of the operators HH, HG, GH, GG, HSH, HSG, GSH, GSG HHS, GHS, HGS, GGS HSHS, HSGS, GSHS and GSGS are applied. Then the same collection of operators is applied to all the derived smooths, i.e. HH, HSH, HHS and HSHS. 

So the next level is obtained from the previous level with basically HH, HG, GH and GG but with extra shifts in the horizontal, vertical and horizontal and vertical directions. The index provides a way to enumerate the paths through this tree where each smooth has 4 children and indexed by a number between 0 and 3. 

Each of the 4 children has 4 components: a smooth, horizontal, vertical and diagonal detail, much in the same way as for the Mallat 2D wavelet transform implemented in the WaveThresh function \code{\link{imwd}}. 
}
\value{
A matrix containing the packet of the 2D non-decimated wavelet coefficients that you require. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{putpacket.wst2D}}, \code{\link{wst2D}}, \code{\link{wst2D.object}}. 
}
\examples{
#
# Create a random image. 
#
myrand <- matrix(rnorm(16), nrow=4, ncol=4)
#myrand
#            [,1]       [,2]        [,3]       [,4] 
#[1,]  0.01692807  0.1400891 -0.38225727  0.3372708
#[2,] -0.79799841 -0.3306080  1.59789958 -1.0606204
#[3,]  0.29151629 -0.2028172 -0.02346776  0.5833292
#[4,] -2.21505532 -0.3591296 -0.39354119  0.6147043
#
# Do the 2D non-decimated wavelet transform
#
myrwst2D <- wst2D(myrand)
#
# Let's access the finest scale detail, not shifted in the vertical
# direction.
#
getpacket(myrwst2D, nlevelsWT(myrwst2D)-1, index=0, type="V")
#           [,1]       [,2] 
#[1,] -0.1626819 -1.3244064
#
# Compare this to the ordinary 2D DWT for the vertical detail at this
# resolution level
imwd(myrand)[[lt.to.name( 1, "DC")]]
#[1] -0.1626819 -1.3244064  1.4113247 -0.7383336
#
# The same numbers but they're not in matrix format because
# imwd returns vectors not matrices.
#
# Now back to the wst2D object. Let's
# extract vertical detail again at level 1 but this time the horizontally
# shifted data.
#
getpacket(myrwst2D, level=1, index=1, type="V")
#           [,1]      [,2] 
#[1,] -0.5984427 0.2599445
#[2,] -0.6502002 1.8027955
#
# So, yes, different data. Now how about at a deeper resolution level.
# Lets have a horizontal shift, as before, for the level 1 but follow it
# with a diagonal shift and this time extract the smooth component:
#
getpacket(myrwst2D, level=0, index=13, type="S")
#           [,1] 
#[1,] -0.5459394
#
# Of course, only one number because this is at level 0
}
\keyword{manip}
\author{G P Nason}
