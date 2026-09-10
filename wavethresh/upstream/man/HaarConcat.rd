\name{HaarConcat}
\alias{HaarConcat}
\usage{
HaarConcat()
}
\arguments{
None
}
\title{Generate a concatenated Haar MA process}
\description{
This function generates a particular set of four concatenated Haar MA
processes.
}
\details{
This function generates a realization of particular kind of non-stationary time series probability model. The returned time series is the result of concatenating 4 time series each of length 128 from the Haar MA process generator
(\code{\link{HaarMA}}) of orders 1, 2, 3 and 4.
The standard deviation of the innovations is 1.
This function was used to generate the figure of the concatenated Haar MA process in Nason, von Sachs and Kroisandt. It produces a kind of time series that can be sparsely represented by the wavelet machinery but at the same time is non-stationary.

See Nason, von Sachs and Kroisandt (2000) Wavelet processes and adaptive
estimation of the evolutionary wavelet spectrum.
\emph{J R Statist Soc, B}, \bold{62}, 271-292.
}
\value{
A vector containing 512 observations from four concatenated Haar MA processes
}
\seealso{\code{\link{HaarMA}},\code{\link{ewspec}}}
\examples{
#
# Generate the concatenated Haar MA process.
#
MyHaarCC <- HaarConcat()
#
# Plot it
#
\dontrun{ts.plot(MyHaarCC)}
}
\author{G P Nason}
\keyword{smooth}
