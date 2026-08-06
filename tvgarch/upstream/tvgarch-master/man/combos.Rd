\name{combos}
\alias{combos}
\alias{combinations}
\title{
  Compute all combinations of a hierarchy of models of n variables, and enumerate the combinations of the elements of a vector.
}

\usage{
  combos(n)
  combinations(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE)
}

\arguments{
  \item{n}{an integer: the number of variables (\code{combos}) or the size of the source vector (\code{combinations})}

  \item{r}{size of the target vectors}

  \item{v}{source vector. Defaults to \code{1:n}}

  \item{set}{logical flag indicating whether duplicates should be removed from
  the source vector \code{v}. Defaults to \code{TRUE}.}

  \item{repeats.allowed}{logical flag indicating whether the constructed
  vectors may include duplicated values. Defaults to \code{FALSE}.}
}

\details{
  \code{combos} lists hierarchy of all possible combinations of n variables in ascending
order, starting with 1 variable, then all combinations of 2 variables,
and so on until the one combination with all n variables. It is used by function
\code{tvgarch} to constrain the size coefficients when s > 1 required to guarantee the variance is positive for all t.

  When using \code{combinations}, the number of combinations increases rapidly with
\code{n} and \code{r}! To use values of \code{n} above about 45, you will need to increase R's
recursion limit. See the \code{expression} argument to the \code{options} command for details
on how to do this. The source code is adapted from the function with the same name
in the package \pkg{gtools}. There, it is stated that the code of the function is
from an email by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help dated Tue,
14 Dec 1999 11:14:04 +0000 (GMT) in response to Alex Ahgarin \email{datamanagement@email.com}.
Original version was named "subsets" and was Written by Bill Venables

}

\value{
  \item{combos:}{a matrix with zeroes in empty elements and 1 in all full elements.}
  \item{combinations:}{a matrix where each row contains a vector of length
\code{r}.}
}

\description{
  \code{combos} produces a matrix of combinations of 1 to n variables in ascending order. \code{combinations} enumerates the possible combinations of a specified size
from the elements of a vector.
}

\examples{

combos(3)

combinations(3,2,letters[1:3])
combinations(3,2,letters[1:3],repeats=TRUE)

}

\references{
Venables, Bill.  "Programmers Note", R-News, Vol 1/1, Jan. 2001.
\url{https://cran.r-project.org/doc/Rnews/}
}

\author{
  \code{combos} by Chris Walsh \email{cwalsh@unimelb.edu.au}, with modifications by Susana 
  Campos-Martins. Original versions of \code{combinations} by Bill Venables
  \email{Bill.Venables@cmis.csiro.au}. Extended to handle \code{repeats.allowed}
  by Gregory R. Warnes \email{greg@warnes.net}.
}

\seealso{
  \code{\link{tvgarch}}
}           
