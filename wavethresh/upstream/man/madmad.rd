\name{madmad}
\alias{madmad}
\title{Compute square of median absolute deviation (mad) function.}
\description{
This function simply returns the square of the median absolute deviation (mad) function in S-Plus. This is required for supply to the \code{\link{threshold}} series of functions which require estimates of spread on the variance scale (not the standard deviation scale). 
}
\usage{
madmad(x)
}
\arguments{
\item{x}{The vector for which you wish to compute the square of mad on.}
}
\value{
The square of the median absolute deviation of the coefficients supplied by \code{x}. 
}
\note{Its a MAD MAD world!}
\section{RELEASE}{Version 3.4.1 Copyright Guy Nason 1994 }
\seealso{
\code{\link{threshold}}
}
\examples{
#
#
# Generate some normal data with mean 0 and sd of 8
# and we'll also contaminate it with an outlier of 1000000
# This is akin to signal wavelet coefficients mixing with the noise.
#
ContamNormalData <- c(1000000, rnorm(1000, mean=0, sd=8))
#
# What is the variance of the data?
#
var(ContamNormalData)
# [1] 999000792
#
# Wow, a seriously unrobust answer!
#
# How about the median absolute deviation?
#
mad(ContamNormalData)
# [1] 8.14832
# 
# A much better answer!
#
# Now let's use madmad to get the answer on the variance scale 
#
madmad(ContamNormalData)
# [1] 66.39512
#
# The true variance was 64, so the 66.39512 was a much better answer
# than that returned by the call to the variance function.
}
\keyword{arith}
\author{G P Nason}
