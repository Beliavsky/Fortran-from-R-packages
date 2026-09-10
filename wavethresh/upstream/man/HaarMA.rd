\name{HaarMA}
\alias{HaarMA}
\title{Generate Haar MA processes.}
\description{
This function generates an arbitrary number of observations from a Haar MA process of any order with a particular variance. 
}
\usage{
HaarMA(n, sd=1, order=5)
}
\arguments{
\item{n}{The number of observations in the realization that you want to create. Note that n does NOT have to be a power of two.}
\item{sd}{The standard deviation of the innovations.} 
\item{order}{The order of the Haar MA process.} 
}
\details{
A Haar MA process is a special kind of time series moving-average (MA) process. A Haar MA process of order k is a MA process of order
\eqn{2^k}. The coefficients of the Haar MA process are given by the filter coefficients of the discrete Haar wavelet at different scales. 

For examples: the Haar MA process of order 1 is an MA process of order 2.
The coefficients are \eqn{1/\sqrt{2}}{1/sqrt(2)} and
\eqn{-1/\sqrt{2}}{-1/sqrt(2)}.
The Haar MA process of order 2 is an MA process of order 4. The coefficients are 1/2, 1/2, -1/2, -1/2 and so on. It is possible to define other processes for other wavelets as well. 

Any Haar MA process is a good examples of a (stationary) LSW process because it is sparsely representable by the locally-stationary wavelet machinery defined in Nason, von Sachs and Kroisandt.

}
\value{
A vector containing a realization of a Haar MA process of the specified order, standard deviation and number of observations. 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern.}
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{HaarConcat}}, \code{\link{ewspec}}, 
}
\examples{
#
# Generate a Haar MA process of order 1 (high frequency series)
#
MyHaarMA <- HaarMA(n=151, sd=2, order=1)
#
# Plot it
#
\dontrun{ts.plot(MyHaarMA)}
#
# Generate another Haar MA process of order 3 (lower frequency), but of
# smaller variance
#
MyHaarMA2 <- HaarMA(n=151, sd=1, order=3)
#
# Plot it
#
\dontrun{ts.plot(MyHaarMA2)}
#
# Let's plot them next to each other so that you can really see the
# differences.
# 
# Plot a vertical dotted line which indicates where the processes are
# joined
#
\dontrun{ts.plot(c(MyHaarMA, MyHaarMA2))}
\dontrun{abline(v=152, lty=2)}
}
\keyword{manip}
\author{G P Nason}
