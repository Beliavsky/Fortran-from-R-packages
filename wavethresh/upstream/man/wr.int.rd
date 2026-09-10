\name{wr.int}
\alias{wr.int}
\title{Computes inverse "wavelets on the interval" transform. }
\description{
This function actually computes the inverse of the "wavelets on the interval" transform. 
}
\usage{
\method{wr}{int}(wav.int.object, \dots)
}
\arguments{
\item{wav.int.object}{A list with components defined by the return from the \code{\link{wd.int}} function.}
\item{\dots}{any other arguments}
}
\details{
(The WaveThresh implementation of the ``wavelets on the interval transform'' was coded by Piotr Fryzlewicz, Department of Mathematics, Wroclaw University of Technology, Poland; this code was largely based on code written by Markus Monnerjahn, RHRK, Universitat Kaiserslautern; integration into WaveThresh by \code{GPN}).
 
See the help on the "wavelets on the interval code" in the wd help page.
}
\value{
The inverse wavelet transform of the wav.int.object supplied. 
}
\note{
It is not recommended that the casual user call this function. The "wavelets on the interval" transform is best called in \code{WaveThresh} via the \code{\link{wd}} function with the argument bc argument set to "\code{interval}". 
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Piotr Fryzlewicz and Markus Monnerjahn 1995-9).}
\seealso{
\code{\link{wd.int}}, \code{\link{wd}}, \code{\link{wr}}. 
}
\examples{
#
# The user is expected to call the wr
# for inverting a "wavelets on the interval transform".
#
}
\keyword{smooth}
\keyword{nonlinear}

\author{Piotr Fryzlewicz and Markus Monnerjahn}
