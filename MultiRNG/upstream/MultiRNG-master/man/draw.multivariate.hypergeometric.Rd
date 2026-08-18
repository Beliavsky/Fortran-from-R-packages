\name{draw.multivariate.hypergeometric}
\alias{draw.multivariate.hypergeometric}
\title{Pseudo-Random Number Generation under Multivariate Hypergeometric Distribution}
\description{This function implements pseudo-random number generation for a multivariate hypergeometric distribution.
}
\usage{
draw.multivariate.hypergeometric(no.row,d,mean.vec,k)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{mean.vec}{Number of items in each category.}
  \item{k}{Number of items to be sampled. Must be a positive integer.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\references{
Demirtas, H. (2004). Pseudo-random number generation in R for commonly used multivariate distributions. \emph{Journal of Modern Applied Statistical Methods}, \bold{3(2)}, 485-497.
}
\examples{
meanvec=c(10,10,12) ; myk=5
mydata=draw.multivariate.hypergeometric(no.row=1e5,d=3,mean.vec=meanvec,k=myk)
apply(mydata,2,mean)-myk*meanvec/sum(meanvec)
}
