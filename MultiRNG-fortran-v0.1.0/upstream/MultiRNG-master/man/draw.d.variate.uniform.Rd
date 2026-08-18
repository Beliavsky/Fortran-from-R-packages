\name{draw.d.variate.uniform}
\alias{draw.d.variate.uniform}
\title{Pseudo-Random Number Generation under Multivariate Uniform Distribution}
\description{This function implements pseudo-random number generation for a multivariate uniform distribution with specified mean vector and covariance matrix.
}
\usage{
draw.d.variate.uniform(no.row,d,cov.mat)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{cov.mat}{Variance-covariance matrix.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\references{
Falk, M. (1999). A simple approach to the generation of uniformly distributed random variables with prescribed correlations. \emph{Communications in Statistics, Simulation and Computation}, \bold{28(3)}, 785-791.
}
\examples{
cmat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
mydata=draw.d.variate.uniform(no.row=1e5,d=3,cov.mat=cmat)
apply(mydata,2,mean)-rep(0.5,3)
cor(mydata)-cmat
}
