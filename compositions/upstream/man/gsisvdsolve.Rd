\name{gsi.svdsolve}
\alias{gsi.svdsolve}
\title{Internal function: Solves singular and non square equation systems}
\description{
  Based on the singular value decomposition, a singular equation system ax=b
  is solved.
}
\usage{
gsi.svdsolve(a,b,...,cond=1E-10)
}
\arguments{
  \item{a}{the matrix of ax=b (a.k.a. left-hand side matrix)}
  \item{b}{the vector or matrix b of ax=b (a.k.a right-hand side, 
             independent element)}
  \item{cond}{the smallest-acceptable condition of the matrix. Smaller
    singular values are truncate}
  \item{\dots}{additional arguments to svd}
}
\value{
  The "smallest" vector or matrix solving this system with minimal 
  joint error among all vectors. 
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
#A <- matrix(c(0,1,0,0,0,0),ncol=2)
#b <- diag(3)
#erg <- gsi.svdsolve(A,b)
#erg
#A \%*\% erg 
#diag(c(0,1,0))  # richtig
}
\keyword{internal}
