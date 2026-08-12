\name{get.partitions}
%\Rdversion{1.0}
\alias{get.partitions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Additive partitioning of natural numbers}
\description{This function solves the problem of additive partitioning of positive integers. The approach for additive partitioning is based on a generating function discussed in details in Voinov and Nikulin (1995). The function enumerates all partitions of a positive integer \code{n} on at most (or exactly) \code{M} parts, \code{M <= n}. 
}
\usage{
 get.partitions(n, M, at.most=TRUE)

}
\arguments{
   \item{n}{A positive integer to be partitioned.}
   \item{M}{A positive integer, the number of parts of \code{n}, \code{M <= n}.}
   \item{at.most}{If \code{TRUE} then partitioning of \code{n} into at most \code{M} parts, if \code{FALSE} then partitioning on exactly \code{M} parts.} 
}

\value{
  \item{p.n}{total number of partitions obtained.}
  \item{partitions}{a matrix with each column presenting partitions of \code{n}.}
}

\author{
Vassilly Voinov, Natalya Pya Arnqvist, Yevgeniy Voinov
}
\references{
  Voinov, V. and Nikulin, M. (1995) Generating functions, problems of additive number theory, and some statistical applications. Revue Roumaine de \enc{Mathématiques}{Mathematiques} Pures et \enc{Appliquées}{Appliquees}, 40(2), 107-147

 Voinov, V.G. and Pya, N.E. (2017) R-software for additive partitioning of positive integers. Mathematical Journal (ISSN 1682-0525) 17(1), 69-76.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{nilde-package}}, \code{\link{get.knapsack}}, \code{\link{get.subsetsum}}, \code{\link{nlde}}     
}

\examples{
## getting all partitions of n = 8 on at most 6 parts...
get.partitions(8,6,at.most=TRUE)

## getting all partitions of n = 8 on exactly 6 parts...
b <- get.partitions(8,6,at.most=FALSE)
b
colSums(b$partitions)

}
\keyword{optimize}% at least one, from doc/KEYWORDS

