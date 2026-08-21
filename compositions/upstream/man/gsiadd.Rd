\name{gsi.add}
\alias{gsi.add}
\alias{gsi.sub}
\alias{gsi.mul}
\alias{gsi.div}
\title{Internal functions: Parallel operations of single and multiple datasets}
\description{
  The given operations are performed in parallel for multiple datasets
  or for two single datasets, or for multiple datasets with a single dataset.

}
\usage{
gsi.add(x,y)
gsi.sub(x,y)
gsi.mul(x,y)
gsi.div(x,y)
          }
\arguments{
  \item{x}{a vector or a matrix}
  \item{y}{a vector or a matrix}
}
\value{
  a vector or a matrix with the operated values 
}
\details{
  All operations +,-,*,/ are performed on unclassed objects.
}

\note{It is better not to use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{perturbe}}
}
\examples{
tmp1 <- matrix(1:12,ncol=3)
tmp2 <- 1:3
#gsi.add(tmp1,tmp2)
#gsi.sub(tmp1,tmp2)
#gsi.mul(tmp1,tmp2)
#gsi.div(tmp1,tmp2)

#gsi.add(tmp2,tmp2)
#gsi.sub(tmp2,tmp2)
#gsi.mul(tmp2,tmp2)
#gsi.div(tmp2,tmp2)

#gsi.add(tmp1,tmp1)
#gsi.sub(tmp1,tmp1)
#gsi.mul(tmp1,tmp1)
#gsi.div(tmp1,tmp1)

}
\keyword{internal}
