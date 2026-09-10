\name{rfftwt}
\alias{rfftwt}
\title{Weight a Fourier series sequence by a set of weights}
\usage{
rfftwt(xrfft, wt)
}
\arguments{
\item{xrfft}{The Fourier series sequence to weight}
\item{wt}{The weights}
}
\description{
Weight the real Fourier series \code{xrfft} of even length by a weight
sequence \code{wt}. The first term of \code{xrfft} is left alone, and the
weights are then applied to pairs of terms in \code{xrfft}. Note:
\code{wt} is half the length of \code{xrfft}.
}
\details{
Description says all
}
\value{
The weighted sequence
}
\seealso{\code{\link{rfft}}}
\author{Bernard Silverman}
\keyword{math}
