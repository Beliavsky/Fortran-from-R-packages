\name{ScalingFunction}
\alias{ScalingFunction}
\title{Compute scaling functions on internally predefined grid}
\usage{
ScalingFunction(filter.number = 10, family = "DaubLeAsymm", resolution = 4096,
    itlevels = 50)
}
\arguments{
\item{filter.number}{The filter number of the associated wavelet. See
	\code{\link{filter.select}}}
\item{family}{The family of the associated wavelet. See
	\code{\link{filter.select}}}
\item{resolution}{The nominal resolution, the actual grid size might be
	larger than this}
\item{itlevels}{The number of complete filtering operations to generate
	the answer}
}
\description{
This is a subsidiary routine not intended to be called by a user:
use \code{\link{draw}} instead.
Generates scaling functions by inserting a Kronecker delta function
into the bottom of the inverse DWT and repeating the inverting steps.
}
\details{
Description says all
}
\value{
A list containing the \code{x} and \code{y} values of the required
scaling function.
}
\seealso{\code{\link{draw}}}
\author{G P Nason}
\keyword{dplot}
