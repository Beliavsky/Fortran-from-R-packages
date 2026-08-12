\name{occurs}
\alias{occurs}
\alias{count}
\title{
  Finding Subsequences
}
\description{
  Counts items, or finds subsequences of (integer) sequences.
}
\usage{
  count(x, sorted = TRUE)

  occurs(subseq, series)
}
\arguments{
  \item{x}{array of items, i.e. numbers or characters.}
  \item{sorted}{logical; default is to sort items beforehand.}

  \item{subseq}{vector of integers.}
  \item{series}{vector of integers.}
}
\details{
  \code{count} counts the items, similar to \code{table}, but as fast and
  a more tractable output. If \code{sorted} then the total number per item
  will be counted, else per repetition.

  If \code{m} and \code{n} are the lengths of \code{s} and \code{S} resp.,
  then \code{occurs(s, S)} determines all positions \code{i} such that
  \code{s == S[i, ..., i+m-1]}.

  The code is vectorized and relatively fast. It is intended to complement
  this with an implementation of Rabin-Karp, and possibly Knuth-Morris-Pratt
  and Boyer-Moore algorithms.
}
\value{
  \code{count} returns a list with components \code{v} the items and
  \code{e} the number of times it apears in the array.
  \code{occurs} returns a vector of indices, the positions where the
  subsequence appears in the series.
}
\examples{
##  Examples
patrn <- c(1,2,3,4)
exmpl <- c(3,3,4,2,3,1,2,3,4,8,8,23,1,2,3,4,4,34,4,3,2,1,1,2,3,4)
occurs(patrn, exmpl)
## [1]  6 13 23

\dontrun{
set.seed(2437)
p <- sample(1:20, 1000000, replace=TRUE)
system.time(i <- occurs(c(1,2,3,4,5), p))  #=>  [1] 799536
##  user  system elapsed 
## 0.017   0.000   0.017 [sec]

system.time(c <- count(p))
##  user  system elapsed 
## 0.075   0.000   0.076 
print(c)
## $v
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
## $e
##  [1] 49904 50216 49913 50154 49967 50045 49747 49883 49851 49893
## [11] 50193 50024 49946 49828 50319 50279 50019 49990 49839 49990
}
}
