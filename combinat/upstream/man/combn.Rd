\name{combn}
\alias{combn}
\alias{combn2}
\alias{nCm}
\title{ Generate all combinations of the elements of x taken m at a time. }
\description{
             Generate all combinations of the elements of x taken m at a time.
             If x is a positive integer,  returns all combinations
             of the elements of seq(x) taken m at a time.
             If argument "fun" is not null,  applies a function given
             by the argument to each point.  If simplify is FALSE,  returns
             a list; else returns a vector or an array.  "..." are passed
             unchanged to function given by argument fun,  if any.

         combn2:Generate all combinations of the elements of x taken two at a time.
         If x is missing,  generate all combinations of 1:n taken two
         at a time (that is,  the indices of x that would give all
         combinations of the elements of x if x with length n had been given).
         Exactly one of arguments "x" and "n" should be given; no provisions
         for function evaluation.

nCm: Compute the binomial coefficient ("n choose m"),  where n is any
real number and m is any integer.  Arguments n and m may be vectors;
they will be replicated as necessary to have the same length.
Argument tol controls rounding of results to integers.  If the
difference between a value and its nearest integer is less than tol,
the value returned will be rounded to its nearest integer.  To turn
off rounding, use tol = 0.  Values of tol greater than the default
should be used only with great caution, unless you are certain only
integer values should be returned.

}
\usage{
combn(x, m, fun=NULL, simplify=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ vector source for combinations }
  \item{m}{ number of elements }
  \item{fun}{ function to be applied to each combination (may be null) }
  \item{simplify}{ logical, if FALSE, returns a list,
otherwise returns vector or array }
  \item{\dots}{ args to fun }
}
\details{
Nijenhuis, A. and Wilf, H.S. (1978) Combinatorial Algorithms for
Computers and Calculators.  NY:  Academic Press.

}
\value{
see simplify argument
}
\references{ ~put references to the literature/web site here ~ }
\author{ Code by Scott Chasalow, R package and doc prep by Vince Carey, stvjc@channing.harvard.edu}
\examples{
 combn(letters[1:4], 2)
 combn(10, 5, min)  # minimum value in each combination
# Different way of encoding points:
 combn(c(1,1,1,1,2,2,2,3,3,4), 3, tabulate, nbins = 4)
#Compute support points and (scaled) probabilities for a
#Multivariate-Hypergeometric(n = 3, N = c(4,3,2,1)) p.f.:
# table.mat(t(combn(c(1,1,1,1,2,2,2,3,3,4), 3, tabulate,nbins=4)))
}
%\keyword{ combinatorics }
\keyword{ models }


