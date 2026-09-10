\name{first.last}
\alias{first.last}
\title{Build a first/last database for wavelet transforms.}
\description{
This function is not intended for user use, but is used by various functions involved in computing and displaying wavelet transforms. It basically constructs "bookeeping" vectors that \code{WaveThresh} uses for working out where coefficient vectors begin and end. 
}
\usage{
first.last(LengthH, DataLength, type, bc="periodic", current.scale=0)
}
\arguments{
\item{LengthH}{Length of the filter used to produce a wavelet decomposition.} 
\item{DataLength}{Length of the data before transforming. This must be a power of 2, say \eqn{2^m}.} 
\item{type}{The type of wavelet transform. Can be "wavelet" or "periodic"}
\item{bc}{This character string argument determines how the boundaries of the the function are to be handled. The permitted values are \code{periodic} or \code{symmetric}. }
\item{current.scale}{Can handle a different initial scale, but usually
	left at the default}
}
\details{
Suppose you begin with \eqn{2^m=2048} coefficients. At the next level you would expect 1024 smoothed data coefficients, and 1024 wavelet coefficients, and if \code{bc="periodic"} this is indeed what happens. However, if \code{bc="symmetric"} you actually need more than 1024 (as the wavelets extend over the edges). The first last database keeps track of where all these "extras" appear and also where they are located in the packed vectors C and D of pyramidal coefficients within wavelet structures. 

For examples, given a \code{first.last.c row} of 

\deqn{-2 3 20}{-2 3 20}

The actual coefficients would be 
\deqn{c_{-2}, c_{-1}, c_{0}, c_{1}, c_{2}, c_{3}}{c_{-2}, c_{-1}, c_{0}, c_{1}, c_{2}, c_{3}}

In other words, there are 6 coefficients, starting at -2 and ending at 3, and the first of these (\eqn{c_{-2}}) appears at an offset of 20 from the beginning of the \code{$C} component vector of the wavelet structure. 

You can ``do'' \code{first.last} in your head for \code{periodic} boundary handling but for more general boundary treatments (e.g. \code{symmetric}) \code{first.last} is indispensable.

}
\value{
A first/last database structure, a list containing the following information: 

\item{first.last.c}{A (m+1)x3 matrix. The first column specifies the real index of the first coefficient of the smoothed data at a level, the 2nd column is the real index of the last coefficient, the last column specifies the offset of the first smoothed datum at that level. The offset is used by the C code to work out where the beginning of the sequence is within a packed vector of the pyramid structure. The first and 2nd columns can be used to work out how many numbers there are at a level. If \code{bc="periodic"} then the pyramid is a true power of 2 pyramid, that is it starts with a power of 2, and the next level is half of the previous. If \code{bc="symmetric"} then the pyramid is nearly exactly a power of 2, but not quite, see the Details section for why this is so. }
\item{ntotal}{The total number of smoothed data/original data points.} 
\item{first.last.d}{A mx3 matrix. As for \code{first.last.c} but for the wavelet coefficients packed as the D component of a wavelet structure.} 
\item{ntotal.d}{The total number of wavelet coefficients.} 
}
\references{Nason, G.P. and Silverman, B.W. (1994). The discrete wavelet transform in S.}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\seealso{
\code{\link{wd}}, \code{\link{wr}}, \code{\link{wr.wd}}, \code{\link{accessC}}, \code{\link{accessD}}, \code{\link{filter.select}}. \code{\link{imwd}}. 
}
\examples{
#
#If you're twisted then you may just want to look at one of these.
#
first.last(length(filter.select(2)), 64)
#$first.last.c:
#First Last Offset
#[1,]     0    0    126
#[2,]     0    1    124
#[3,]     0    3    120
#[4,]     0    7    112
#[5,]     0   15     96
#[6,]     0   31     64
#[7,]     0   63      0
#
#$ntotal:
#[1] 127
#
#$first.last.d:
#First Last Offset
#[1,]     0    0     62
#[2,]     0    1     60
#[3,]     0    3     56
#[4,]     0    7     48
#[5,]     0   15     32
#[6,]     0   31      0
#
#$ntotal.d:
#[1] 63
#
#
}
\keyword{manip}
\author{G P Nason}
