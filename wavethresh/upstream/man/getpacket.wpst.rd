\name{getpacket.wpst}
\alias{getpacket.wpst}
\title{Get packet of coefficients from a non-decimated wavelet packet object (wpst). }
\description{
This function extracts and returns a packet of coefficients from a non-decimated wavelet packet (\code{\link{wpst}}) object. 
}
\usage{
\method{getpacket}{wpst}(wpst, level, index, \dots )
}
\arguments{
\item{wpst}{Non-decimated wavelet packet object from which you wish to extract the packet from.} 
\item{level}{The resolution level of the coefficients that you wish to extract. Can range from 0 to \code{\link{nlevelsWT}}(wpst). The coefficients at level \code{\link{nlevels}} are the data the created the \code{wpst.object}. 
}
\item{index}{The index number within the resolution level of the packet of coefficients that you wish to extract. Index ranges from 0 to
\eqn{(4^r)-1} where \code{r = nlevelsWT - level}. }
\item{\dots}{any other arguments}
}
\details{
The \code{\link{wpst}} transform produces a non-decimated wavelet packet object. This is a "cross" between a \code{wavelet packet} object and a \code{non-decimated wavelet} object. In other words the transform produces \emph{wavelet packet} coefficients at every possible integer shift (unlike the ordinary wavelet packet transform which is aligned to a dyadic grid). 

Each packet of coefficients is obtained by chaining together the effect of the two \emph{packet operators} DG and DH: these are the high and low pass quadrature mirror filters of the Mallat pyramid algorithm scheme followed by both even \emph{and} odd decimation. For a full description of this algorithm and how coefficients are stored within see Nason, Sapatinas and Sawczenko, 1998. 

Note that this function extracts \emph{packets}. If you want to obtain the wavelet packet coefficients for each shift you need to use the \code{\link{accessD.wpst}}function. This function extracts particular wavelet packet coefficients for a particular shift. In particular, this function returns a number of coefficients dependent on the scale level requested whereas \code{\link{accessD.wpst}} always returns a vector of coefficients of length equal to the input data that created the \code{wpst.object} initially. 
}
\value{
A vector containing the packet of non-decimated wavelet packet coefficients that you wished to extract. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{accessD.wpst}}, \code{\link{wpst}}, 
}
\examples{
#
# Create some random data
#
myrand <- rnorm(16)
#myrand
# [1]  0.19268626 -0.41737181 -0.30806613  0.07435407  0.99871757
# [6] -0.58935121 -1.38049759 -0.13346631  1.55555403 -1.60581265
#[11]  0.14353621  1.21277774  1.13762337 -1.08577934 -0.29745609
#[16]  0.50977512
#
# Do the non-decimated wavelet packet transform
#
myrwpst <- wpst(myrand)
#
# Let's access what is a level nlevelsWT(myrwpst)
#
getpacket(myrwpst, nlevelsWT(myrwpst), index=0)
# [1]  0.19268626 -0.41737181 -0.30806613  0.07435407  0.99871757
# [6] -0.58935121 -1.38049759 -0.13346631  1.55555403 -1.60581265
#[11]  0.14353621  1.21277774  1.13762337 -1.08577934 -0.29745609
#[16]  0.50977512
#
# I.e. the data that created the object.
#
# How about extracting the 3rd (last) packet at level 3?
#
getpacket(myrwpst, 3, index=3)
#[1] -2.660657144  0.688415755 -1.764060698  0.717267105 -0.206916242
#[6] -0.659983747  0.005836952 -0.196874007
#
# Of course, there are only 8 coefficients at this level.
}
\keyword{manip}
\author{G P Nason}
