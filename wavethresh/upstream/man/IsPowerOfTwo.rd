\name{IsPowerOfTwo}
\alias{IsPowerOfTwo}
\title{Decides whether vector elements are integral powers of two (returns NA if not).
}
\description{
This function checks to see whether its input is a power of two. If it is then it returns that power otherwise it returns NA. 
}
\usage{
IsPowerOfTwo(n)
}
\arguments{
\item{n}{Vector of numbers that are to be checked whether it is a power of two.}
}
\details{
Function takes the log of the input, divides this by log(2) and if the result is integral then it knows the input is true power of two. 
}
\value{
If \code{n} is a power of two, then the power is returned otherwise \code{NA} is returned. 
}
\section{RELEASE}{Version 3.6.0 Copyright Guy Nason 1995 }
\seealso{
\code{\link{nlevelsWT.default}}. 
}
\examples{
#
# Try and see whether 1,2,3 or 4 are powers of two! 
#
IsPowerOfTwo(1:4)
# [1]  0  1 NA  2
#
# Yes, 1,2 and 4 are the 0, 1 and 2nd power of 2. However, 3 is not an
# integral power of two.
}
\keyword{arith}
\author{G P Nason}
