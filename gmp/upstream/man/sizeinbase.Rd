\name{sizeinbase}
\alias{sizeinbase}

\title{Compute size of a bigz in a base}
\description{
 Return an approximation to the number of
 character the integer X would have printed in base b.  The
 approximation is never too small. 

 In case of powers of 2, function gives exact result.
}
\usage{
sizeinbase(a, b=10)
}

\arguments{
  \item{a}{big integer, i.e. \code{"\link{bigz}"}}
  \item{b}{base}
}

\value{
  integer of the same length as \code{a}: the size, i.e. number of
  digits, of each \code{a[i]}.
}
\references{The GNU MP Library, see \url{https://gmplib.org}}
\author{Antoine Lucas}

\examples{
sizeinbase(342434, 10)# 6 obviously

Iv <- as.bigz(2:7)^500
sizeinbase(Iv)
stopifnot(sizeinbase(Iv)       == nchar(as.character(Iv)),
          sizeinbase(Iv, b=16) == nchar(as.character(Iv, b=16)))
}
\keyword{arith}

