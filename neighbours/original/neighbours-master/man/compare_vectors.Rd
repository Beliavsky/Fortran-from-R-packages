\name{compare_vectors}
\alias{compare_vectors}
\title{
  Compare Vectors
}
\description{
  Compare numeric or logical vectors.
}
\usage{
compare_vectors(..., sep = "", diff.char = "|")
}
\arguments{
  \item{\dots}{
    vectors of the same length
  }
  \item{sep}{
    a string
  }
  \item{diff.char}{
    a single character
  }
}
\details{

  The function compares vectors with one another. The
  main purpose is to print a useful representation of
  differences (and return differences, usually
  invisibly).

  The function is still experimental and will likely change.

}
\value{

  depends on how the function is called; typically a list

}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
  Enrico Schumann
}
\seealso{
  \code{\link{neighbourfun}}
}
\examples{
x <- runif(5) > 0.5
nb <- neighbourfun(type = "logical")

compare_vectors(x, nb(x))
## 01010
##  |
## 00010
## The vectors differ in  1  place.

nb <- neighbourfun(type = "logical", stepsize = 2)
compare_vectors(x, nb(x))
## 01010
## |   |
## 11011
## The vectors differ in  2  places.
}
