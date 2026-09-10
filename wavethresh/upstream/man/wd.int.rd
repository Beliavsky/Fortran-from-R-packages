\name{wd.int}
\alias{wd.int}
\title{Computes "wavelets on the interval" transform}
\description{
This function actually computes the "wavelets on the interval" transform. 

\bold{NOTE:} It is not recommended that the casual user call this function. The "wavelets on the interval" transform is best called in \code{WaveThresh} via the \code{\link{wd}} function with the argument bc argument set to \code{"interval"}.
 }
\usage{
wd.int(data, preferred.filter.number, min.scale, precond)
}
\arguments{
\item{data}{The data that you wish to apply the "wavelets on the interval" transform to.}
\item{preferred.filter.number}{Which wavelet to use to do the transform. This is an integer ranging from 1 to 8. See the Cohen, Daubeches and Vial (1993) paper. Wavelets that do not "overlap" a boundary are just like the ordinary Daubechies' wavelets.}
\item{min.scale}{At which resolution level to transform to.}
\item{precond}{If true performs preconditioning of the input vector to try and ensure that simple polynomial sequences (less than in order to the wavelet used) map to zero elements.}
}
\details{
(The \code{WaveThresh} implementation of the ``wavelets on the interval transform'' was coded by Piotr Fryzlewicz, Department of Mathematics, Wroclaw University of Technology, Poland; this code was largely based on code written by Markus Monnerjahn, RHRK, Universitat Kaiserslautern; integration into WaveThresh by GPN). 

See the help on the "wavelets on the interval code" in the \code{\link{wd}} help page. 
}
\value{
A list containing the wavelet transform of the \code{data}. We again emphasize that this list is not intended for human consumption, use the \code{\link{wd}} function with the correct \code{bc="interval"} argument. 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Piotr Fryzlewicz and Markus Monnerjahn 1995-9). }
\seealso{
\code{\link{wd}}, \code{\link{wr}}, \code{\link{wr.int}}. 
}
\examples{
#
# The user is expected to call the wr
# for inverting a "wavelets on the interval transform" and not to use
# this function explicitly
#
}
\keyword{smooth}
\keyword{nonlinear}
\author{Piotr Fryzlewicz}
