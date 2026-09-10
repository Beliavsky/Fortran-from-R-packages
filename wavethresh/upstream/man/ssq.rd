\name{ssq}
\alias{ssq}
\title{Compute sum of squares difference between two vectors}
\usage{
ssq(u,v)
}
\arguments{
\item{u}{One of the vectors}
\item{v}{The other of the vectors}
}
\description{
Given two vectors, u and v, of length n, this function computes
\eqn{\sum_{i=1}^n (u_i - v_i)^2}{sum((u-v)^2)}.
}
\details{
Description says all
}
\value{
The sum of squares difference between the two vectors
}
\examples{
ssq(c(1,2), c(3,4))
#[1] 8
}
\author{G P Nason}
\keyword{math}
