\name{mfirst.last}
\alias{mfirst.last}
\title{Build a first/last database for multiple wavelet transforms. }
\description{
This function is not intended for user use, but is used by various functions involved in computing and displaying multiple wavelet transforms.
}
\usage{
mfirst.last(LengthH, nlevels, ndecim, type = "wavelet", bc = "periodic")
}
\arguments{
\item{LengthH}{Number of filter matrix coefficients.} 
\item{nlevels}{Number of levels in the decomposition} 
\item{ndecim}{The decimation scale factor for the multiple wavelet basis.} 
\item{type}{Whether the transform is non-decimated or ordinary (wavelet). The non-decimated multiple wavelet transform is not yet supported.}
\item{bc}{This argument determines how the boundaries of the the function are to be handled. The permitted values are periodic or \code{symmetric} }
}
\details{
Suppose you begin with \code{2^m}=2048 coefficient vectors. At the next level you would expect 1024 smoothed data vectors, and 1024 wavelet vectors, and if \code{bc="periodic"} this is indeed what happens. However, if \code{bc="symmetric"} you actually need more than 1024 (as the wavelets extend over the edges). The first last database keeps track of where all these "extras" appear and also where they are located in the packed vectors C and D of pyramidal coefficients within wavelet structures. 

For examples, given a \code{first.last.c} row of 

\deqn{-2 3 20}{-2 3 20} 

The `position' of the coefficient vectors would be 

\deqn{c_{-2}, c_{-1}, c_{0}, c_{1}, c_{2}, c_{3}}{c_{-2}, c_{-1}, c_{0}, c_{1}, c_{2}, c_{3}}

In other words, there are 6 coefficients, starting at -2 and ending at 3, and the first of these (\eqn{c_{-2}}) appears at column 20 of the \code{$C} component matrix of the wavelet structure. 

You can ``do'' first.last in your head for periodic boundary handling but for more general boundary treatments (e.g. symmetric) first.last is indispensable. 

The numbers in first last databases were worked out from inequalities derived from: Daubechies, I. (1988). 
}
\value{
A first/last database structure, a list containing the following information: 

\item{first.last.c}{A \code{(m+1)x3} matrix. The first column specifies the real index of the first coefficient vector of the smoothed data at a level, the 2nd column is the real index of the last coefficient vector, the last column specifies the offset of the first smoothed datum at that level. The offset is used by the C code to work out where the beginning of the sequence is within a packed vector of the pyramid structure. The first and 2nd columns can be used to work out how many numbers there are at a level. If bc="periodic" then the pyramid is a true power of 2 pyramid, that is it starts with a power of 2, and the next level is half of the previous. If bc="symmetric" then the pyramid is nearly exactly a power of 2, but not quite, see the Details section for why this is so.}
\item{nvecs.c}{The number of C coefficient vectors.} 
\item{first.last.d}{A \code{mx3} matrix. As for \code{first.last.c} but for the wavelet coefficients packed as the D component of a wavelet structure.}
\item{nvecs.d}{The number of \code{D} coefficient vectors.} 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6)}
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mwd.object}}, \code{\link{mwd}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{threshold.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
}
\examples{
#
#To see the housekeeping variables for a decomposition with
# 4 filter coefficient matices
# 5 resolution levels and a decimation scale of two
# use:
mfirst.last(4,5,2)
# $first.last.c:
#	 First Last Offset
# [1,]     0    0     62
# [2,]     0    1     60
# [3,]     0    3     56
# [4,]     0    7     48
# [5,]     0   15     32
# [6,]     0   31      0
# 
# $nvecs.c:
# [1] 63
# 
# $first.last.d:
#	 First Last Offset
# [1,]     0    0     30
# [2,]     0    1     28
# [3,]     0    3     24
# [4,]     0    7     16
# [5,]     0   15      0
# 
# $nvecs.d:
# [1] 31
}
\keyword{datagen}
\author{Tim Downie}
