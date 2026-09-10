\name{AvBasis.wst}
\alias{AvBasis.wst}
\title{Perform basis averaging for (packet-ordered) non-decimated wavelet transform.}
\description{
Perform basis averaging for (packet-ordered) non-decimated wavelet transform.
}
\usage{
\method{AvBasis}{wst}(wst, Ccode=TRUE, \dots)
}
\arguments{
\item{wst}{An object of class \code{\link{wst}} that contains coefficients of a packet ordered non-decimated wavelet transform (e.g. produced by the \code{\link{wst}} function.}
\item{Ccode}{If TRUE then fast compiled C code is used to perform the transform. If FALSE then S code is used. Almost always use the default TRUE option. (It is conceivable that some implementation can not use the C code and so this option permits use of the slower S code).}
\item{\dots}{any other arguments}
}
\details{
The packet-ordered non-decimated wavelet transform computed by \code{\link{wst}} computes the coefficients of an input vector with respect to a library of all shifts of wavelet basis functions at all scales. Here "all shifts" means all integral shifts with respect to the finest scale coefficients, and "all scales" means all dyadic scales from 0 (the coarsest) to J-1 (the finest) where \code{2^J = n} where \code{n} is the number of data points of the input vector. As such the packet-ordered non-decimated wavelet transform contains a library of all possible shifted wavelet bases.
 
\code{Basis selection} It is possible to select a particular basis and invert that particular representation. In WaveThresh a basis is selected by creating a \code{nv} (node.vector) class object which identifies the basis. Then the function \code{\link{InvBasis}} takes the wavelet representation and the node.vector and inverts the representation with respect to the selected basis. The two functions \code{\link{MaNoVe}} and \code{\link{numtonv}} create a node.vector: the first by using a \code{Coifman-Wickerhauser} minimum entropy best-basis algorithm and the second by basis index. 

\bold{Basis averaging}. Rather than select a basis it is often useful to preserve information from all of the bases. For examples, in curve estimation, after \link{threshold}ing a wavelet representation the coefficients are coefficients of an estimate of the truth with respect to all of the shifted basis functions. Rather than select \emph{one} of them we can average over all estimates. This sometimes gives a better curve estimate and can, for examples, get rid of Gibbs effects. See Coifman and Donoho (1995) for more information on how to do curve estimation using the packet ordered non-decimated wavelet transform, thresholding and basis averaging. 

Further it might seem that inverting each wavelet transform and averaging might be a computationally expensive operation: since each wavelet inversion costs order \eqn{n} operations and there are n different bases and so you might think that the overall order is \eqn{n^2}.
It turns out that since many of the coarser scale basis functions are duplicated between bases there is redundancy in the non-decimated transform. Coifman and Donoho's TI-denoising algorithm makes use of this redundancy which results in an algorithm which only takes order \eqn{n \log n}{n*logn} operations. 

For an examples of denoising using the packet-ordered non-decimated wavelet transform and basis averaging see Johnstone and Silverman, 1997. The WaveThresh implementation of the basis averaging algorithm is to be found in Nason and Silverman, 1995 
}
\value{
A vector containing the average of the wavelet representation over all the basis functions. The length of the vector is \code{2^nlev} where \code{nlev} is the number of levels in the input \code{wst} object.
}
\section{RELEASE}{Version 3.6.0 Copyright Guy Nason 1995}
\seealso{
\code{\link{av.basis}},
\code{\link{wst}}, \code{\link{wst.object}}, \code{\link{MaNoVe}},  \code{\link{numtonv}}, \code{\link{InvBasis}}, \code{\link{wavegrow}}
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
#
# Now take the packet-ordered non-decimated wavelet transform
#
tdwst <- wst(test.data)
#
# Now "invert" it using basis averaging
#
tdwstAB <- AvBasis(tdwst)
#
# Let's compare it to the original
#
sum( (tdwstAB - test.data)^2)
#
# [1] 9.819351e-17 
#
# Very small. They're essentially same.
#
# See the threshold.wst help page for an
# an examples of using basis averaging in curve estimation.
}
\keyword{manip} 
\author{G P Nason}
