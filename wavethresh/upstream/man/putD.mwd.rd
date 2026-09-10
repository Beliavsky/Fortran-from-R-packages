\name{putD.mwd}
\alias{putD.mwd}
\title{Put wavelet coefficients into multiple wavelet structure }
\description{
The wavelet coefficients from a multiple wavelet decomposition structure, \code{\link{mwd.object}}, (e.g. returned from \code{\link{mwd}}) are packed into a single matrix in that structure. This function copies the \code{\link{mwd.object}}, replaces some wavelet coefficients in the copy, and then returns the copy. 
}
\usage{
\method{putD}{mwd}(mwd, level, M, boundary = FALSE, index = FALSE, \dots)
}
\arguments{
\item{mwd}{Multiple wavelet decomposition structure whose coefficients you wish to replace.}
\item{level}{The level that you wish to replace.}
\item{M}{Matrix of replacement coefficients.}
\item{boundary}{If \code{boundary} is \code{FALSE} then only the "real" data is replaced (and it is easy to predict the required length of \code{M}). If \code{boundary} is \code{TRUE} then you can replace the boundary values at a particular level as well (but it is hard to predict the required length of\code{M}, and the information has to be obtained from the \code{mfirst.last} database component of \code{mwd}).}
\item{index}{If index is \code{TRUE} then the index numbers into the \code{mwd$D} array where the matrix \code{M} would be stored is returned. Otherwise, (default) the modified \code{\link{mwd.object}} is returned. }
\item{\dots}{any other arguments}
}
\details{
The \code{\link{mwd}} function produces a wavelet decomposition structure. 

The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as a linear matrix of coefficients. PutD obtains information about where the wavelet coefficients appear from the fl.dbase component of mwd, in particular the array \code{fl.dbase$first.last.d} which gives a complete specification of index numbers and offsets for \code{mwd$D}. 

Note also that this function only puts information into mwd class objects. To extract coefficients from mwd structures you have to use the accessD.mwd function. 

See Downie and Silverman, 1998. 
}
\value{
An object of class \code{\link{mwd.object}} if index is \code{FALSE}, otherwise the index numbers indicating where the \code{M} matrix would have been inserted into the \code{mwd$D} object are returned. 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6). }
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Generate an mwd object
#
tmp <- mwd(rnorm(32))
#
# Now let's examine the finest resolution detail...
#
accessD(tmp, level=3)
#          [,1]       [,2]       [,3]      [,4]       [,5]        [,6] 
#[1,] 0.8465672  0.4983564  0.3408087 0.1340325  0.5917774 -0.06804291
#[2,] 0.6699962 -0.2535760 -1.0344445 0.2068644 -0.4912086  1.16039885
#           [,7]       [,8] 
#[1,] -0.6226445  0.2617596
#[2,] -0.4956576 -0.5555795
#
#
# A matrix. There are two rows one for each mother wavelet in this 
# two-ple multiple wavelet transform and at level 3 there are 2^3 columns.
#
# Let's set the coefficients of the first mother wavelet all equal to zero
# for this examples
#
newdmat <- accessD(tmp, level=3)
newdmat[1,] <- 0
#
# Ok, let's insert it back at level 3
#
tmp2 <- putD(tmp, level=3, M=newdmat)
#
# And check it
#
accessD(tmp2, level=3)   
#          [,1]      [,2]      [,3]      [,4]       [,5]     [,6]       [,7] 
#[1,] 0.0000000  0.000000  0.000000 0.0000000  0.0000000 0.000000  0.0000000
#[2,] 0.6699962 -0.253576 -1.034445 0.2068644 -0.4912086 1.160399 -0.4956576
#           [,8] 
#[1,]  0.0000000
#[2,] -0.5555795
#
#
# Yep, all the first mother wavelet coefficients at level 3 are now zero.
}
\keyword{manip}
\author{Tim Downie}
