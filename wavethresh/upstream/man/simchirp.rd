\name{simchirp}
\alias{simchirp}
\title{Compute and return simulated chirp function.}
\description{
This function computes and returns the coordinates of the reflected simulated chirp function described in Nason and Silverman, 1995. This function is a useful test function for evaluating wavelet shrinkage and time-scale analysis methodology as its frequency changes over time.
}
\usage{
simchirp(n=1024)
}
\arguments{
\item{n}{The number of ordinates from which to sample the chirp signal.} 
}
\details{
This function computes and returns the x and y coordinates of the reflected chirp function described in Nason and Silverman, 1995. 

The formula for the reflected simulated chirp is  *formula*

The chirp returned is a discrete sample on \code{n} equally spaced points between -1 and 1. 

}
\value{
A list with two components: 

\item{x}{a vector of length \code{n} containing the ordered x ordinates of the chirp from -1 to 1.}
\item{y}{a vector of length \code{n} containing the corresponding y ordinates of the chirp.} 
}
\section{RELEASE}{Version 3.5.3 Copyright Guy Nason 1994 }
\examples{
#
# Generate the chirp
#
test.data <- simchirp()$y
\dontrun{ts.plot(test.data)}
}
\keyword{manip}
\author{G P Nason}
