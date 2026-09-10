\name{permn}
\alias{permn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generates all permutations of the elements of x }
\description{
Generates all permutations of the elements of x, in a minimal-
change order. If x is a positive integer,  returns all permutations
of the elements of seq(x). If argument "fun" is not null,  applies
a function given by the argument to each point. "..." are passed
unchanged to the function given by argument fun, if any.
}
\usage{
permn(x, fun=NULL, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ vector }
  \item{fun}{ if non.null, applied at each perm }
  \item{\dots}{ args passed to fun }
}
%\details{
%}
\value{
list: each component is either a permutation, or the
results of applying fun to a permutation
}
\references{ 
Reingold, E.M., Nievergelt, J., Deo, N. (1977) Combinatorial
Algorithms: Theory and Practice. NJ: Prentice-Hall. pg. 170.
}
\seealso{ 
sample, fact, combn, hcube, xsimplex
}

\examples{
    # Convert output to a matrix of dim c(6, 720)
     t(array(unlist(permn(6)), dim = c(6, gamma(7))))
       # A check that every element occurs the same number of times in each
       # position
       apply(t(array(unlist(permn(6)), dim = c(6, gamma(7)))), 2, tabulate,
               nbins = 6)

       # Apply, on the fly, the diff function to every permutation
       t(array(unlist(permn(6, diff)), dim = c(5, gamma(7))))
}
\keyword{ models }

