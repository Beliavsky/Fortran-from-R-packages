\name{Shannon.entropy}
\alias{Shannon.entropy}
\title{Compute Shannon entropy}
\description{
Computes Shannon entropy of the squares of a set of coefficients. 
}
\usage{
Shannon.entropy(v, zilchtol=1e-300)
}
\arguments{
\item{v}{A vector of coefficients (e.g. wavelet coefficients).} 
\item{zilchtol}{A small number. Any number smaller than this is considered to be zero for the purposes of this function.} 
}
\details{
This function computes the Shannon entropy of the squares of a set of coefficients. The squares are used because we are only interested in the entropy of the energy of the coefficients, not their actual sign. 

The entropy of the squares of \code{v} is given by \code{sum( v^2 * log(v^2) )}. In this implementation any zero coefficients (determined by being less than \code{zilchtol}) have a zero contribution to the entropy. 

The Shannon entropy measures how "evenly spread" a set of numbers is. If the size of the entries in a vector is approximately evenly spread then the Shannon entropy is large. If the vector is sparsely populated or the entries are very different then the Shannon entropy is near zero. Note that the input vectors to this function usually have their norm normalized so that diversity of coefficients corresponds to sparsity. 
}
\value{
A number representing the Shannon entropy of the input vector. 
}
\section{RELEASE}{Version 3.7.2 Copyright Guy Nason 1996 }
\seealso{
\code{\link{MaNoVe.wst}}, \code{\link{wst}}, 
}
\examples{
#
# Generate some test data
#
#
# A sparse set
#
Shannon.entropy(c(1,0,0,0))
#0
#
# A evenly spread set
#
Shannon.entropy( rep( 1/ sqrt(4), 4 ))
#1.386294
}
\keyword{manip}
\author{G P Nason}
