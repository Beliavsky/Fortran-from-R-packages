\name{AvBasis.wst2D}
\alias{AvBasis.wst2D}
\title{Perform basis averaging for (packet-ordered) 2D non-decimated wavelet transform.}
\description{
Perform basis averaging for (packet-ordered) 2D non-decimated wavelet transform. 
}
\usage{
\method{AvBasis}{wst2D}(wst2D, \dots)
}
\arguments{
\item{wst2D}{An object of class \code{\link{wst2D}} that contains coefficients of a packet ordered 2D non-decimated wavelet transform (e.g. produced by the \code{\link{wst2D}} function.}
\item{\dots}{any other arguments}
}
\details{
The packet-ordered 2D non-decimated wavelet transform computed by \code{\link{wst2D}} computes the coefficients of an input matrix with respect to a library of all shifts of wavelet basis functions at all scales. Here "all shifts" means all integral shifts with respect to the finest scale coefficients with shifts in both the horizontal and vertical directions, and "all scales" means all dyadic scales from 0 (the coarsest) to J-1 (the finest) where \code{2^J = n} where \code{n} is the dimension of the input matrix. As such the packet-ordered 2D non-decimated wavelet transform contains a library of all possible shifted wavelet bases.
 
\bold{Basis averaging}. Rather than select \emph{a} basis it is often useful to preserve information from all of the bases. For examples, in curve estimation, after thresholding, the coefficients are coefficients of an estimate of the truth with respect to all of the shifted basis functions. Rather than select one of them we can average over all estimates. This sometimes gives a better curve estimate and can, for examples, get rid of Gibbs effects. See Coifman and Donoho (1995) for more information on how to do curve estimation using the packet ordered non-decimated wavelet transform, thresholding and basis averaging. See Lang et al. (1995) for further details of surface/image estimation using the 2D non-decimated DWT. 
}
\value{
A square matrix of dimension $2^nlevelsWT$ containing the average-basis ``reconstruction'' of the \code{\link{wst2D}} object. 
}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998}
\seealso{
\code{\link{wst2D}}, \code{\link{wst2D.object}}
}
\examples{
#
# Generate some test data
#
#test.data <- matrix(rnorm(16), 4,4)
#
# Now take the 2D packet ordered DWT 
#
#tdwst2D <- wst2D(test.data)
#
# Now "invert" it using basis averaging
#
#tdwstAB <- AvBasis(tdwst2D)
#
# Let's compare it to the original
#
#sum( (tdwstAB - test.data)^2)
#
# [1] 1.61215e-17
#
# Very small. They're essentially same.
#
}
\keyword{manip}
\author{G P Nason}
