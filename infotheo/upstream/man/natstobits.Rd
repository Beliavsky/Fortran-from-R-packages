\name{natstobits}
\alias{natstobits}

\title{convert nats into bits}
\usage{natstobits(H)}
\arguments{
  \item{H}{double denoting a value (in nats), as returned by one of the function of the infotheo package}
}
\value{ \code{natstobits} returns a double that is the conversion of the nats value into bits.}
\description{
	  \code{natstobits} takes a value in nats (a double) as input and returns the value in bits (a double). 
}
\details{
Information-theoretic quantities can have different units depending on the base of the logarithm used in their computation.
All the function of tha package use a base e, hence the unit is the nat. The value in bit is given by using the base 2, 
hence the conversion is done by multiplying by log2(e) = 1.442695.
}
\author{
  Patrick E. Meyer
}

\examples{
  data(USArrests)
  H <- entropy(discretize(USArrests))
  natstobits(H)
}
\keyword{misc}
