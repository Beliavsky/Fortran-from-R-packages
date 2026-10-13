\name{slideWeight}
\alias{slideWeight}
\title{ Generate Time Weights Flexiibly}
\description{
Returns a vector of weights between 1 and 0 (inclusive).
The general form is: some recent weights are 1, some early weights are 0, 
and weights between these vary linearly.
It is allowable for no points to be specified as 0 or 1.
}
\usage{
slideWeight(n, fractions = c(0, 1), observations = NULL, 
	locations = NULL)
}
\arguments{
  \item{n}{ 
	the length of the result -- the number of observations.
}
  \item{fractions}{ 
	length 2 numeric vector giving the fraction of the way 
	between 0 and \code{n} where the ends of the linear slide 
	should be.
	There is no requirement that the fractions need to be between 0 and 1.

	This is ignored if \code{observations} or \code{locations} is given.
 }
  \item{observations}{ 
	length 2 numeric vector giving the number of observations that 
	should be 1 and the number of observations that should have 
	any weight.
	The smaller is taken to be the former and the larger as the latter.
 }
  \item{locations}{ 
	length 2 numeric vector giving the locations (subscripts) for the
	ends of the linear slide.
	The smaller number is the index of the last 0, 
	the larger is the index of the last value smaller than 1.
 }
}
\value{
	a numeric vector of length \code{n} with values within the range
	of 0 and 1 (inclusive).
	The values are non-decreasing -- meaning put more weight on more
	recent observations.
}

\section{Details }{
The arguments \code{fractions}, \code{observations} and \code{locations}
are all to control where the ends of the weights strictly between 0 and 1 are.
Only one of them is used.
\code{locations} takes precedence.
\code{fractions} is used if neither of the other two are given.

This is suitable to give as the \code{weights} argument to 
\code{var.shrink.eqcor} and \code{factor.model.stat}.
}

\section{Revision }{
This help was last revised 2014 March 09.
}

\seealso{ 
\code{\link{var.shrink.eqcor}}, \code{\link{factor.model.stat}}.
}
\examples{
# examples assume number of observations is 200

# all weights either 0 or 1
rollwin50 <- slideWeight(200, observations=c(50,50))

# 50 with full weight, 100 more with partial weight
swght50.150 <- slideWeight(200, observations=c(50, 150))

# approximately default weights of var.shrink.eqcor and factor.model.stat
lindec3 <- slideWeight(200, fractions=c(-1/2, 3/2))
}
\keyword{ ts }
