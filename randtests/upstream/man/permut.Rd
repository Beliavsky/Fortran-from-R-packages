\name{permut}
\alias{permut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate all permutations of \eqn{m} elements of a vector
}
\description{
Generate all permutations of \code{x} taken \eqn{m} at a time. If argument FUN is not NULL, applies a function given by the argument to each permutation. 
}
\usage{
permut(x, m=length(x), FUN=NULL,...)
}
\arguments{
  \item{x}{vector source for permutations.}
  \item{m}{number of elements to choose. Default is \code{m=length(x)}.}
  \item{FUN}{function to be applied to each permutation; default NULL means the identity, i.e., to return the permutation.}
  \item{...}{optionally, further arguments to FUN.}
}
\details{
Based on function \code{permutations} from package \code{gtools}. This function is required for the computation of the exact p-value of some randomness tests.
}
\value{
A matrix with one permutation, or the value returned by FUN, in each line.
}
%\references{
%}
%\examples{
%}
%\keyword{  }
%\keyword{ randomness test }% __ONLY ONE__ keyword per line
