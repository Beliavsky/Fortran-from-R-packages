\name{zlag}
\alias{zlag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the lag of a vector.}
\description{Computes the lag of a vector, with missing elements replaced
by NA}
\usage{
zlag(x, d= 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector}
  \item{d}{compute the lag d of x}
}
\value{
A vector whose k-th element equals x[k-d] with x[t]=NA for t<=0}
\author{Kung-Sik Chan}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
x=1:5
zlag(x,2)
}
\keyword{methods}
