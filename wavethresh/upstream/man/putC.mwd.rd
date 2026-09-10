\name{putC.mwd}
\alias{putC.mwd}
\title{Put smoothed data into wavelet structure }
\description{
The smoothed and original data from a multiple wavelet decomposition structure, \code{\link{mwd.object}}, (e.g. returned from \code{\link{mwd}}) are packed into a single matrix in that structure. This function copies the \code{\link{mwd.object}}, replaces some smoothed data in the copy, and then returns the copy.
}
\usage{
\method{putC}{mwd}(mwd, level, M, boundary = FALSE, index = FALSE, \dots)
}
\arguments{
\item{mwd}{Multiple wavelet decomposition structure whose coefficients you wish to replace.}
\item{level}{The level that you wish to replace.}
\item{M}{Matrix of replacement coefficients.}
\item{boundary}{If \code{boundary} is \code{FALSE} then only the "real" data is replaced (and it is easy to predict the required length of \code{M}). If boundary is TRUE then you can replace the boundary values at a particular level as well (but it is hard to predict the required length of \code{M}, and the information has to be obtained from the \code{mfirst.last} database component of \code{mwd}).}
\item{index}{If index is \code{TRUE} then the index numbers into the \code{mwd$C} array where the matrix \code{M} would be stored is returned. Otherwise, (default) the modified \code{\link{mwd.object}} is returned.}
\item{\dots}{any other arguments}
}
\details{
The \code{\link{mwd}} function produces a wavelet decomposition structure. 

The need for this function is a consequence of the pyramidal structure of Mallat's algorithm and the memory efficiency gain achieved by storing the pyramid as a linear matrix of coefficients. PutC obtains information about where the smoothed data appears from the fl.dbase component of mwd, in particular the array \code{fl.dbase$first.last.c} which gives a complete specification of index numbers and offsets for \code{mwd$C}. 

Note also that this function only \emph{puts} information into \code{\link{mwd}} class objects. To \emph{extract} coefficients from \code{\link{mwd}} structures you have to use the \code{\link{accessC.mwd}} function. 

See Downie and Silverman, 1998. 
}
\value{
An object of class \code{\link{mwd.object}} if \code{index} is \code{FALSE}, otherwise the index numbers indicating where the \code{M} matrix would have been inserted into the \code{mwd$C} object are returned. 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6). }
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
# Generate an mwd object
#
tmp <- mwd(rnorm(32))
#
# Now let's examine the finest resolution smooth...
#
accessC(tmp, level=3)
#           [,1]       [,2]       [,3]       [,4]        [,5]      [,6] 
#[1,] -0.4669103 -1.3150580 -0.7094966 -0.1979214  0.32079986 0.5052254
#[2,] -0.7645379 -0.8680941  0.1004062  0.6633268 -0.05860848 0.5757286
#          [,7]       [,8] 
#[1,] 0.5187380  0.6533843
#[2,] 0.2864293 -0.4433788
#
# A matrix. There are two rows one for each father wavelet in this 
# two-ple multiple wavelet transform and at level 3 there are 2^3 columns.
#
# Let's set the coefficients of the first father wavelet all equal to zero
# for this examples
#
newcmat <- accessC(tmp, level=3)
newcmat[1,] <- 0
#
# Ok, let's insert it back at level 3
#
tmp2 <- putC(tmp, level=3, M=newcmat)
#
# And check it
#
accessC(tmp2, level=3)   
#           [,1]       [,2]      [,3]      [,4]        [,5]      [,6]      [,7] 
#[1,]  0.0000000  0.0000000 0.0000000 0.0000000  0.00000000 0.0000000 0.0000000
#[2,] -0.7645379 -0.8680941 0.1004062 0.6633268 -0.05860848 0.5757286 0.2864293
#           [,8] 
#[1,]  0.0000000
#[2,] -0.4433788
#
# Yep, all the first father wavelet coefficients at level 3 are now zero.
}
\keyword{manip}
\author{G P Nason}
